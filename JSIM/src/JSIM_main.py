'''Code for the bulk framework of the JWST/NIRSpec
simulator. This code should show the main functions/processes
to move from an input datacube (lambda, y, x) to output cubes:
 - Observed IFU datacube corresponding to mock observation
 - Background datacube (or spectrum) containing all BG sources
 - Noise cube

Written by Simon Zieleniewski

Started 16-10-15

Last edited 03-12-15
'''

#Import all required modules
import numpy as n
import astropy.io.fits as p
import scipy.constants as sc
import os
import pprint as pp
import multiprocessing as mp
from JSIM_Specres import spectral_res
from JSIM_Lowres import low_res_mode
from JSIM_Loop import pp_wavelength_loop, wavelength_loop
from JSIM_Background import create_background_cube
from JSIM_Transmission import create_thruput_cube
from modules.generate_cubes import *
from modules.spaxel_rebin import *
from modules.misc_utils import generate_seed, point_source
from modules.fits_utils import *
from JSIM_Config import *



def main(datacube, outdir, DIT, NDIT, grating, ignoreLSF, res_jitter=0,
         site_temp=280.5, combine_ndits=True, Spec_nyquist=True, Spec_samp=1.,
         noise_force_seed=0, remove_background='False', return_object='False',
         return_transmission='False', drizzle=0, version='1', nprocs=mp.cpu_count()-1):
    '''Main body of code: gathers all inputs and calls all functions
    to generate outputs.

    Inputs:
        datacube: Input high resolution datacube (RA, DEC, lambda)
        outdir: Output file directory
        DIT: Exposure time [s]
        NDIT: No. of exposures
        grating: Spectral grating
        ignoreLSF: turn off LSF convolution in spectral dimension
        res_jitter: Residual telescope jitter [mas]
        site_temp: Telescope temperature [K]
        combine_ndits: Boolean - Combine NDITS
        spec_Nyquist: Boolean - Set spectral sampling to Nyquist sample spectral resolution element
        Spec_samp: Spectral sampling [A/pix] - used if spec_Nyquist = False
        noise_force_seed: Force random number seed to take a set value (0=No, 1-10=yes and takes that value)
        remove_background: Boolean - Subtract background spectrum
        return_object: Boolean - Return object cube
        return_transmission: Boolean - Return throughput cube
        drizzle: Boolean - Output sampling of 50 mas
        version: Version of code. Equal to SVN revision number'
        nprocs: no. of processes to use for parallel processing

    Outputs:
        observed_cube: FITS file of datacube representing HARMONI observation of input cube
        background_cube: FITS file of datacube representing background sources
        reduced_cube: FITS file of observed datacube after subtracting background
        noise_cube: FITS file of noise for each pixel in datacube

    '''

    print 'Filename: ', datacube
    print 'Output dir: ', outdir
    print 'DIT = ', DIT
    print 'NINT = ', NDIT
    print 'Grating = ', grating
    print 'Ignore LSF?', ignoreLSF
    print 'Residual jitter = ', res_jitter
    print 'Temperature = ', site_temp
    print 'combine NINTS? ', combine_ndits
    print 'Noise seed = ', noise_force_seed
    print 'Spectral Nyquist sampling? ', Spec_nyquist
    print 'Subtract background? ', remove_background
    print 'Return Object cube?', return_object
    print 'Return transmission?', return_transmission
    print 'Drizzle?', drizzle
    print 'No. of processes = ', nprocs

    #Generate random number seed
    seednum = generate_seed(noise_force_seed)
    n.random.seed(seednum)

    #Spaxels
    if drizzle:
        spax = (50., 50.)
    else:
        spax = config_data['spaxels']


    #OPEN INPUT FITS file
    if os.path.isfile(datacube) == True and os.path.splitext(datacube)[1].lower() == '.fits':
        cube, head = p.getdata(datacube, 0, header=True, memmap=True)
    #If not FITS file, try passing datacube directly
    else:
        try:
            cube = datacube.data
            head = datacube.header
        except:
            raise ValueError('UNKNOWN INPUT FILE FORMAT: Please use FITS file')


    #Check that datacube has required headers to be processed in simulator
    #Required headers = ['CDELT1/2/3'], ['CRVAL3'], ['NAXIS1/2/3'], ['FUNITS'], ['CRPIX3'] = 1,
    #['CTYPE1/2/3'] = 'RA, DEC, WAVELENGTH, ['CUNIT1/2/3'] = MAS, MAS, microns/angstroms/etc,
    #['SPECRES'].
    fits_header_check(head)
    pp.pprint(head)


    #CREATE WAVELENGTH ARRAY IN MICRONS
    lambs, head = wavelength_array(head)


    #RESCALE DATACUBE TO CHOSEN SPECTRAL RESOLUTION.
    if grating == 'Prism':
        print 'Prism low resolution mode chosen'
        cube, head, lambs, delta_lambda, pix_disp, R = low_res_mode(cube, head, lambs)
    else:
        cube, head, lambs, delta_lambda, pix_disp, R = spectral_res(cube, head, grating, lambs, config_data['gratings'],
                                                                    spec_nyquist=Spec_nyquist, spec_samp=Spec_samp,
                                                                    ignoreLSF=ignoreLSF)


    #POINT SOURCE CHECK
    if (head['NAXIS1']*head['CDELT1'])/float(spax[0]) < 1. or (head['NAXIS2']*head['CDELT2'])/float(spax[1]) < 1.:
        print 'Point source - single spaxel output'
        config_data['ADR'] = 'OFF'
        cube, head = point_source(cube, head, spax)


    #REBIN CUBE TO 10 MAS
    print 'Input datacube sampling = (%.1f, %.1f) mas' % (head['CDELT1'], head['CDELT2'])
    print 'Input PSF sampling = (%.1f, %.1f) mas' % (spax[0]/10.,spax[1]/10.)
    if head['CDELT1'] != spax[0] and head['CDELT2'] != spax[1]:
        print 'Rebinning input datacube to %.0f mas (same as PSF generation scale)' % (spax[0]/10.)
        cube *= (head['CDELT1']*head['CDELT2']*1.E-6)
        cube, head = spaxel_scale(cube, head, (spax[0]/10.,spax[0]/10.))
        cube /= (head['CDELT1']*head['CDELT2']*1.E-6)


    #Empty output cube
    out_cube = n.zeros((len(lambs),int(head['CDELT2']*head['NAXIS2']/spax[1]),
                        int(head['CDELT1']*head['NAXIS1']/spax[0])), dtype=n.float64)

    print 'Input spaxel scale (x, y)', (head['CDELT1'], head['CDELT2']), ' mas'
    print 'Cube shape (lam, y, x) = ', cube.shape
    print 'Output spaxel scale (x, y) = ', spax, ' mas'
    print 'Output cube shape (lam, y, x) = ', out_cube.shape


    #WAVELENGTH CHANNEL LOOP
    print 'Entering loop over wavelength channels'
    ncpus = nprocs
    print 'No. of CPUs: ', ncpus
    try:
        import pprocess
        print 'Using ', str(ncpus), ' CPUs'
    except:
        print 'No pprocess module'
        print 'Using 1 CPU'
        ncpus = 1
    if ncpus >= 3:
        out_cube, head = pp_wavelength_loop(cube, head, lambs, out_cube,
                                            ((head['NAXIS1']*head['CDELT1'])/float(spax[0]),head['NAXIS2']*head['CDELT2']/float(spax[1])),
                                            spax, usecpus=ncpus)
    elif ncpus < 3:
        print 'Using 1 CPU'
        out_cube, head = wavelength_loop(cube, head, lambs, out_cube,
                                         ((head['NAXIS1']*head['CDELT1'])/float(spax[0]),head['NAXIS2']*head['CDELT2']/float(spax[1])),
                                         spax)
    else:
        raise ValueError('Something went wrong at Loop stage!')

    inspax = n.array([head['CDELT1'],head['CDELT2']])*1.E-3
    outspax = n.array([spax[0],spax[1]])*1.E-3
    head['CDELT1'] = (spax[0], "mas")
    head['CDELT2'] = (spax[1], "mas")
    head['NAXIS1'] = int(head['NAXIS1']*head['CDELT1']/float(spax[0]))
    head['NAXIS2'] = int(head['NAXIS2']*head['CDELT2']/float(spax[1]))

    #Remove input cube from memory!
    cube = n.zeros((1,1,1))
    print 'PSF convolution, ADR effect, spaxel rebinning - all done!'


    #Energy-to-Photons Conversion factor will depend on head['FUNITS'] value
    if head['FUNITS'] == 'J/s/m2/um/arcsec2':
        print 'Flux units = ', head['FUNITS']
        en2ph_conv_fac = (sc.h * sc.c)/(lambs*1.E-6) #J
    elif head['FUNITS'] == 'erg/s/cm2/A/arcsec2':
        print 'Flux units = ', head['FUNITS']
        en2ph_conv_fac = (sc.h * sc.c * 1.E7)/(lambs*1.E-6 * 1.E4 * 1.E4) #erg/1.E4(cm2->m2)/1.E4(A->um)
    elif head['FUNITS'] == 'J/s/m2/A/arcsec2':
        print 'Flux units = ', head['FUNITS']
        en2ph_conv_fac = (sc.h * sc.c)/(lambs*1.E-6 * 1.E4) #J/1.E4(A->um)
    elif head['FUNITS'] == 'erg/s/cm2/um/arcsec2':
        print 'Flux units = ', head['FUNITS']
        en2ph_conv_fac = (sc.h * sc.c * 1.E7)/(lambs*1.E-6 * 1.E4) #erg/1.E4(cm2->m2)
    else:
        raise ValueError('UNKNOWN FLUX UNITS: Please change FUNITS header key to erg/s/cm2/A/arcsec2')
    en2ph_conv_fac.shape = (len(lambs),1,1)


    #Total throughput cube
    throughput_cube = create_thruput_cube(out_cube.shape, lambs, delta_lambda, grating, telescope=True, instrument=True)

    #Create object cube (enter spaxel as arcsec = mas*1.E-3)
    #[units of passed values: lambs: [um], DIT: [s], outspax: [arcsec], pix_disp: [um/pixel], area: [m2]]
    out_cube = n.divide(out_cube, en2ph_conv_fac)
    object_cube, object_shot_noise, noiseless_object = generate_object_cube(out_cube, lambs, throughput_cube, DIT,
                                                          NDIT, outspax, pix_disp, config_data['area'])

    #Instrument throughput cube
    inst_qe_cube = create_thruput_cube(out_cube.shape, lambs, delta_lambda, grating, telescope=False, instrument=True)

    #Create background cube [sky + telescope + instrument photons]
    background_cube, background_shot_noise, noiseless_background = create_background_cube(out_cube.shape, lambs, throughput_cube, inst_qe_cube,
                                                                 DIT, NDIT, outspax, delta_lambda, pix_disp, config_data['area'], site_temp,
                                                                 sky=True, telescope=False, emis=config_data['emissivity'])
    #Create observed cube [source + background + dark + SQRT(NDIT)*read]
    read_cube, nrc = generate_read_cube(out_cube.shape, lambs, DIT, NDIT, readsig_vis=config_data['read_noise'],
                                        readsig_nirs=[config_data['read_noise'],config_data['read_noise']])
    dark_cube, noiseless_dark = generate_dark_cube(out_cube.shape, lambs, DIT, NDIT,
                                                   vis_dark=config_data['dark_current'], nir_dark=config_data['dark_current'])
    observed_cube = object_cube + background_cube + dark_cube + read_cube
    observed_noise = generate_noise_cube(noiseless_object, noiseless_background, noiseless_dark, nrc)

    #Create separate background measurement cube
    background_cube2, background_shot_noise2, noiseless_background2 = create_background_cube(out_cube.shape, lambs, throughput_cube, inst_qe_cube,
                                                                   DIT, NDIT, outspax, delta_lambda, pix_disp, config_data['area'], site_temp,
                                                                   sky=True, telescope=False, emis=config_data['emissivity'])
    dark_cube2, noiseless_dark2 = generate_dark_cube(out_cube.shape, lambs, DIT, NDIT,
                                                   vis_dark=config_data['dark_current'], nir_dark=config_data['dark_current'])
    read_cube2, nrc2 = generate_read_cube(out_cube.shape, lambs, DIT, NDIT, readsig_vis=config_data['read_noise'],
                                          readsig_nirs=[config_data['read_noise'],config_data['read_noise']])
    background_cube2 = background_cube2 + dark_cube2 + read_cube2
    background_noise2 = generate_noise_cube(n.zeros_like(noiseless_object), noiseless_background2, noiseless_dark2, nrc2)


    #Return cubes: Observed and background each with noise, or reduced (observed-background) with noise.
    #Return as FITS files

    #Common headers to add to all cubes
    common_header_info = {'DIT': DIT, 'NDIT': NDIT, 'R': R, 'Grating': grating,
                          'RES_JITT': res_jitter, 'FUNITS': 'electrons', 'VERSION': version}

    outFile_objcube = datacube.split('/')[-1].split('.fits')[0] + '_Object_cube'
    outFile_nlobjcube = datacube.split('/')[-1].split('.fits')[0] + '_Noiseless_Object_cube'
    outFile_obscube = datacube.split('/')[-1].split('.fits')[0] + '_Observed_cube'
    outFile_bgrcube = datacube.split('/')[-1].split('.fits')[0] + '_Background_cube'
    outFile_nbgrcube = datacube.split('/')[-1].split('.fits')[0] + '_Noiseless_Background_cube'
    outFile_snr = datacube.split('/')[-1].split('.fits')[0] + '_Obs_SNR_cube'
    outFile_redcube = datacube.split('/')[-1].split('.fits')[0] + '_Reduced_cube'
    outFile_redsnr = datacube.split('/')[-1].split('.fits')[0] + '_Red_SNR_cube'
    outFile_trans = datacube.split('/')[-1].split('.fits')[0] + '_Transmission_cube'


    #return observed_cube, observed_noise, background_cube2, background_noise2
    generate_fits_cube(observed_cube, head, lambs, outFile_obscube, common_header_info, outdir, varext=n.power(observed_noise,2))
    generate_fits_cube(background_cube2, head, lambs, outFile_bgrcube, common_header_info, outdir, varext=n.power(background_noise2,2))
    generate_fits_cube((noiseless_background2+noiseless_dark2), head, lambs, outFile_nbgrcube, common_header_info, outdir)
    generate_fits_cube((noiseless_object/observed_noise), head, lambs, outFile_snr, common_header_info, outdir)

    if return_object == 'True':
        generate_fits_cube(object_cube, head, lambs, outFile_objcube, common_header_info, outdir)
        generate_fits_cube(noiseless_object, head, lambs, outFile_nlobjcube, common_header_info, outdir)

    if remove_background == 'True':
        reduced_cube = n.subtract(observed_cube, background_cube2)
        reduced_noise = generate_noise_cube(noiseless_object, 2.*noiseless_background,
                                            2.*noiseless_dark, n.sqrt(2.)*nrc2)

        generate_fits_cube(reduced_cube, head, lambs, outFile_redcube, common_header_info, outdir, varext=n.power(reduced_noise,2))
        generate_fits_cube((noiseless_object/reduced_noise), head, lambs, outFile_redsnr, common_header_info, outdir)

    if return_transmission == 'True':
        generate_fits_cube(throughput_cube, head, lambs, outFile_trans, common_header_info, outdir)

    print 'Output cubes created!'
    print '          --------     '
    print 'Inbuilt Telescope Parameters'
    print '          --------     '
    print 'Telescope area = %.2f [m^2]' % config_data['area']
    print 'Telescope emissivity = %.2f ' % config_data['emissivity']
    print 'Telescope temperature = %.1f [K]' % site_temp
    print '          --------     '
    print 'Inbuilt Instrument Parameters'
    print '          --------     '
    print 'Wavelength range = %.4g [A] to %.4g [A]' % (lambs[0]*10000., lambs[-1]*10000.)
    print 'Spectral resolving power = %.0i' % R
    if R == 100.:
        print 'Spectral resolution = variable'
        print 'Spectral sampling = 2 pixels per FWHM'
    else:
        print 'Spectral resolution = %.3f [A]' % (delta_lambda*10000.)
        print 'Spectral sampling = %.3f [A/pix]' % (pix_disp*10000.)
    print 'Spatial scale = (%.1f mas, %.1f mas)' % (spax[0], spax[1])
    print 'Detector dark current = %.4f [e-/s/pixel]' % config_data['dark_current']
    print 'Detector read noise = %.4f [e-/pixel] RMS [DIT>120s]' % config_data['read_noise']
    print '                     = %.4f [e-/pixel] RMS [DIT<=120s]' % config_data['read_noise']
    print ' '
    print 'Simulation Complete!'
