'''Code for the bulk framework of the HARMONI
simulator. This code should show the main functions/processes
to move from an input datacube (lambda, y, x) to output cubes:
 - Observed IFU datacube corresponding to mock observation
 - Background datacube (or spectrum) containing all BG sources
 - Noise cube

Written by Simon Zieleniewski

Started 28-05-13
Last edited 27-10-16 by Laurence Routledge
'''

#Import all required modules
import numpy as n
import astropy.io.fits as p
import scipy.constants as sc
import os
import pprint as pp
import multiprocessing as mp
from TheSimulator_Specres import spectral_res
from TheSimulator_Lowres import low_res_mode
from TheSimulator_PSFs import psf_setup
from TheSimulator_ADR import calc_ref, calc_adr, optimalguide
from TheSimulator_Loop import pp_wavelength_loop, wavelength_loop
from TheSimulator_Background import create_background_cube
from TheSimulator_Transmission import create_thruput_cube
from modules.generate_cubes import *
from modules.misc_utils import generate_seed, point_source
from modules.fits_utils import *
###Config file containing hardwired constants:
# - mirror diameter
# - mirror obscuration ratio
# - telescope area
# - telescope emissivity
# - dark current
# - read noise
# - transmission: set to inbuilt or specific constant value
# - ADR ON or OFF
from TheSimulator_Config import *



def main(datacube, outdir, DIT, NDIT, grating, spax, seeing, zenith_ang, telescope='E-ELT', user_PSF='None',
         AO='SCAO', res_jitter=0, site_temp=280.5, combine_ndits=True, Spec_nyquist=True, Spec_samp=1.,
         noise_force_seed=0, remove_background='False', return_object='False', return_transmission='False',
         adr_switch='True', version='1', nprocs=mp.cpu_count()-1):
    '''Main body of code: gathers all inputs and calls all functions
    to generate outputs.

    Inputs:
        datacube: Input high resolution datacube (RA, DEC, lambda)
        outdir: Output file directory
        DIT: Exposure time [s]
        NDIT: No. of exposures
        grating: Spectral grating
        spax: spatial pixel (spaxel) scale [mas]
        seeing: Atmospheric seeing FWHM [arcsec]
        zenith_ang: Zenith angle [deg]
        telescope: Telescope type (E-ELT, VLT)
        user_PSF: Path to user uploaded PSF file
        AO: AO mode (LTAO, SCAO, Gaussian)
        res_jitter: Residual telescope jitter [mas]
        site_temp: Telescope temperature [K]
        combine_ndits: Boolean - Combine NDITS
        spec_Nyquist: Boolean - Set spectral sampling to Nyquist sample spectral resolution element
        Spec_samp: Spectral sampling [A/pix] - used if spec_Nyquist = False
        noise_force_seed: Force random number seed to take a set value (0=No, 1-10=yes and takes that value)
        remove_background: Boolean - Subtract background spectrum
        return_object: Boolean - Return object cube
        return_transmission: Boolean - Return throughput cube
        adr_switch: Boolean - turn ADR on or off.
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
    print 'X spaxels = ', spax[0]
    print 'Y spaxels = ', spax[1]
    print 'Grating = ', grating
    print 'Telescope: ', telescope
    print 'AO = ', AO
    print 'Seeing = ', seeing
    print 'Zenith angle = ', zenith_ang
    print 'Residual jitter = ', res_jitter
    print 'Temperature = ', site_temp
    print 'combine NINTS? ', combine_ndits
    print 'Noise seed = ', noise_force_seed
    print 'Spectral Nyquist sampling? ', Spec_nyquist
    print 'Subtract background? ', remove_background
    print 'User PSF? ' , user_PSF
    print 'Return Object cube?', return_object
    print 'Return transmission?', return_transmission
    print 'ADR off?', adr_switch
    print 'No. of processes = ', nprocs

    #Generate random number seed
    seednum = generate_seed(noise_force_seed)
    n.random.seed(seednum)


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

    #telescope options = 'VLT', 'E-ELT'
    #Sets diameter D, M1 collecting area and obscuration ratio eps
    if telescope == 'VLT':
        D = config_data['VLT']['diam']
        area = config_data['VLT']['area']
        eps = config_data['VLT']['obsc']
    elif telescope == 'E-ELT':
        D = config_data['E-ELT']['diam']
        area = config_data['E-ELT']['area']
        eps = config_data['E-ELT']['obsc']
    else:
        raise ValueError('UNKNOWN TELESCOPE CHOICE: choose E-ELT or VLT')

    #Seeing value at zenith distance calculated as:
    #X = 1/(n.cos(n.radians(zenith_ang)))
    #seeing = seeing*X***(3/5.)
    X = 1./(n.cos(n.radians(zenith_ang)))
    seeing = seeing*X**(3./5.)

    if AO == 'LTAO' or AO == 'GLAO' or AO == 'SCAO':
        if seeing > 1.1:
            print "AO PSFs currently don't support seeing > 1.1'' FWHM"
            seeing = 1.1
        if seeing < 0.67:
            print "AO PSFs currently don't support seeing < 0.67'' FWHM"
            seeing = 0.67
    print 'Seeing at chosen zenith angle = %.3f arcsec FWHM' % seeing


    #CREATE WAVELENGTH ARRAY IN MICRONS
    lambs, head = wavelength_array(head)
    #Check limits of wavelength array and cut cube into HARMONI wavelength range if neccessary
    if grating=='None' and lambs[0] < config_data['gratings']['V+R'][0]:
        print 'Cube lowest wavelength lower than HARMONI range'
        llim = n.where(lambs > config_data['gratings']['V+R'][0])[0][0]
        lambs = lambs[llim:]
        cube = cube[llim:,:,:]
        head['CRVAL3'] = lambs[0]
    if grating=='None' and lambs[-1] > config_data['gratings']['H+K'][1]:
        print 'Cube longest wavelength longer than HARMONI range'
        hlim = n.where(lambs < config_data['gratings']['H+K'][1])[0][-1]
        lambs = lambs[:hlim+1]
        cube = cube[:hlim+1,:,:]
    head['NAXIS3'] = len(lambs)


    #RESCALE DATACUBE TO CHOSEN SPECTRAL RESOLUTION.
    if grating == 'Iz+J+H+K':
        print 'R = 500 low resolution mode chosen'
        cube, head, lambs, delta_lambda, pix_disp, R = low_res_mode(cube, head, lambs)
    else:
        cube, head, lambs, delta_lambda, pix_disp, R = spectral_res(cube, head, grating, lambs, config_data['gratings'],
                                                                 spec_nyquist=Spec_nyquist, spec_samp=Spec_samp)

    #PSF GENERATION CODE
    cube, head, psfspax, psfparams, psfsize, upsf, upsflams = psf_setup(cube, head, lambs, spax,
                                                                        user_PSF, AO, seeing, [D,eps])

    #POINT SOURCE CHECK
    if (head['NAXIS1']*head['CDELT1'])/float(spax[0]) < 1. or (head['NAXIS2']*head['CDELT2'])/float(spax[1]) < 1.:
        print 'Point source - single spaxel output'
        config_data['ADR'] = 'OFF'
        cube, head = point_source(cube, head, spax)

    #CALCULATE ADR
    airmass = 1./n.cos(n.radians(zenith_ang))
    refr = calc_ref(lambs, site_temp)
    optlam = optimalguide(lambs[0], lambs[-1], site_temp)
    adr = calc_adr(lambs, refr, airmass, optlam)
    if adr_switch == 'True':
        print 'Turning ADR off'
        config_data['ADR'] = 'OFF'

    #Empty output cube
    out_cube = n.zeros((len(lambs),int(head['CDELT2']*head['NAXIS2']/spax[1]),
                        int(head['CDELT1']*head['NAXIS1']/spax[0])), dtype=n.float64)


    print 'Cube shape (lam, y, x) = ', cube.shape
    print 'Chosen spaxel scales (x, y) = ', spax
    print 'Output cube shape (lam, y, x) = ', out_cube.shape
    print 'PSF generation size = ', psfsize

    #WAVELENGTH CHANNEL LOOP
    print 'Entering loop over wavelength channels'
    psfvars = [psfparams, psfspax, psfsize, [D, eps], res_jitter, seeing, upsf, upsflams]

    ncpus = nprocs
    print 'No. of CPUs: ', ncpus
    try:
        import pprocess
        print 'Using ', str(ncpus), ' CPUs'
    except:
        print 'No pprocess module'
        print 'Using 1 CPU'
        ncpus = 1
    PSFs = ['LTAO', 'SCAO', 'Gaussian']
    if AO in PSFs and ncpus >= 3:
        out_cube, head = pp_wavelength_loop(cube, head, lambs, out_cube, AO, psfvars, adr,
                                            (head['NAXIS1']*head['CDELT1']/float(spax[0]),head['NAXIS2']*head['CDELT2']/float(spax[1])),
                                            spax, adr_switch=config_data['ADR'], usecpus=ncpus)
    elif AO in PSFs and ncpus < 3:
        # print 'Using 1 CPU'
        out_cube, head = wavelength_loop(cube, head, lambs, out_cube, AO, psfvars, adr,
                                            (head['NAXIS1']*head['CDELT1']/float(spax[0]),head['NAXIS2']*head['CDELT2']/float(spax[1])),
                                            spax, adr_switch=config_data['ADR'])
    # elif AO == 'Gaussian':
    #     out_cube, head = pp_wavelength_loop(cube, head, lambs, out_cube, AO, psfvars, adr,
    #                                         (head['NAXIS1']*head['CDELT1']/float(spax[0]),head['NAXIS2']*head['CDELT2']/float(spax[1])),
    #                                         spax, adr_switch=config_data['ADR'])
    else:
        raise ValueError('UNKNOWN AO CHOICE: choose LTAO, SCAO or Gaussian')

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
    throughput_cube = create_thruput_cube(out_cube.shape, lambs, delta_lambda, grating, zenith_ang, 
    [config_data['trans_w_grat'],config_data['trans_wo_grat']], sky=True, telescope=True, instrument=True, QE=True)

    #Create object cube (enter spaxel as arcsec = mas*1.E-3)
    #[units of passed values: lambs: [um], DIT: [s], outspax: [arcsec], pix_disp: [um/pixel], area: [m2]]
    out_cube = n.divide(out_cube, en2ph_conv_fac)
    object_cube, object_shot_noise, noiseless_object = generate_object_cube(out_cube, lambs, throughput_cube, DIT,
                                                          NDIT, outspax, pix_disp, area)

    #Instrument + Quantum efficiency cube
    inst_qe_cube = create_thruput_cube(out_cube.shape, lambs, delta_lambda, grating,
    [config_data['trans_w_grat'],config_data['trans_wo_grat']], sky=False, telescope=False, instrument=True, QE=True)
    #Quantum efficiency cube
    qe_cube = create_thruput_cube(out_cube.shape, lambs, delta_lambda, grating,
    [config_data['trans_w_grat'],config_data['trans_wo_grat']], sky=False, telescope=False, instrument=False, QE=True)

    #Create background cube [sky + telescope + instrument photons]
    background_cube, background_shot_noise, noiseless_background = create_background_cube(out_cube.shape, lambs, throughput_cube, inst_qe_cube,
                                                                 qe_cube, DIT, NDIT, outspax, delta_lambda, pix_disp, area, site_temp, sky=True,
                                                                 telescope=True, instrument=False, emis=config_data['emissivity'])
    #Create observed cube [source + background + dark + SQRT(NDIT)*read]
    read_cube, nrc = generate_read_cube(out_cube.shape, lambs, DIT, NDIT, readsig_vis=config_data['read_noise_vis'],
                                        readsig_nirs=[config_data['read_noise_nir'],config_data['read_noise_nir_lowexp']])
    dark_cube, noiseless_dark = generate_dark_cube(out_cube.shape, lambs, DIT, NDIT,
                                                   vis_dark=config_data['dark_current_vis'], nir_dark=config_data['dark_current_nir'])
    observed_cube = object_cube + background_cube + dark_cube + read_cube
    observed_noise = generate_noise_cube(noiseless_object, noiseless_background, noiseless_dark, nrc)

    #Create separate background measurement cube
    background_cube2, background_shot_noise2, noiseless_background2 = create_background_cube(out_cube.shape, lambs, throughput_cube, inst_qe_cube,
                                                                   qe_cube, DIT, NDIT, outspax, delta_lambda, pix_disp, area, site_temp, sky=True,
                                                                   telescope=True, instrument=False, emis=config_data['emissivity'])
    dark_cube2, noiseless_dark2 = generate_dark_cube(out_cube.shape, lambs, DIT, NDIT,
                                                   vis_dark=config_data['dark_current_vis'], nir_dark=config_data['dark_current_nir'])
    read_cube2, nrc2 = generate_read_cube(out_cube.shape, lambs, DIT, NDIT, readsig_vis=config_data['read_noise_vis'],
                                          readsig_nirs=[config_data['read_noise_nir'],config_data['read_noise_nir_lowexp']])
    background_cube2 = background_cube2 + dark_cube2 + read_cube2
    background_noise2 = generate_noise_cube(n.zeros_like(noiseless_object), noiseless_background2, noiseless_dark2, nrc2)


    #Return cubes: Observed and background each with noise, or reduced (observed-background) with noise.
    #Return as FITS files

    #Common headers to add to all cubes
    common_header_info = {'DIT': DIT, 'NDIT': NDIT, 'SEEING': seeing, 'ZENITH': zenith_ang, 'AOMODE': AO,
                          'USER_PSF': user_PSF, 'R': R, 'Grating': grating, 'RES_JITT': res_jitter,
                          'FUNITS': 'electrons', 'VERSION': version}

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
    print 'Telescope diameter = %.1f [m]' % D
    print 'Telescope obscuration ratio = %.2f' % eps
    print 'Telescope area = %.2f [m^2]' % area
    print 'Telescope emissivity = %.2f ' % config_data['emissivity']
    print 'Telescope temperature = %.1f [K]' % site_temp
    print 'AO mode = ', AO
    print '          --------     '
    print 'Inbuilt Instrument Parameters'
    print '          --------     '
    print 'Wavelength range = %.4g [A] to %.4g [A]' % (lambs[0]*10000., lambs[-1]*10000.)
    print 'Spectral resolving power = %.0i' % R
    if R == 500.:
        print 'Spectral resolution = variable'
        print 'Spectral sampling = 2 pixels per FWHM'
    else:
        print 'Spectral resolution = %.3f [A]' % (delta_lambda*10000.)
        print 'Spectral sampling = %.3f [A/pix]' % (pix_disp*10000.)
    print 'Spatial scale = (%.1f mas, %.1f mas)' % (spax[0], spax[1])
    print 'Detector dark current = %.4f [e-/s/pixel] for visible' % config_data['dark_current_vis']
    print '                       = %.4f [e-/s/pixel] for near-IR' % config_data['dark_current_nir']
    print 'Detector read noise = %.3f [e-/pixel] RMS for visible' % config_data['read_noise_vis']
    print '                     = %.3f [e-/pixel] RMS for near-IR [DIT>120s]' % config_data['read_noise_nir']
    print '                     = %.3f [e-/pixel] RMS for near-IR [DIT<=120s]' % config_data['read_noise_nir_lowexp']
    print '                      Cut-off at 0.8 um.'
    print ' '
    print 'Simulation Complete!'




if __name__=="__main__":

    import sys

    print 'No. of arguments = ', len(sys.argv)

    if len(sys.argv) != 21:

        print""
        print '    --------     '
        print 'HARMONI Simulator'
        print '    --------     '
        print""
        print 'COMMAND LINE USAGE'
        print ""
        print 'Enter command line arguments in following order:'
        print '1. datacube: Input datacube filepath'
        print '2. DIT: Exposure time [s]'
        print '3. NDIT: No. of exposures'
        print '4. grating: [V+R, Iz+J, H+K, V, R, Iz, J, H, K, V-high, R-high, z, J-high, H-high, K-high, None]'
        print '5. x-spax: x spatial pixel (spaxel) scale [mas]'
        print '6. y-spax: y spatial pixel (spaxel) scale [mas]'
        print '7. seeing: Atmospheric seeing FWHM [arcsec] - Between 0.67"-1.10"'
        print '8. zenith_ang: Zenith angle [deg]'
        print '9. telescope: Telescope type (E-ELT, VLT)'
        print '10. user_PSF: Path to user uploaded PSF file - Enter None if not required'
        print '11. AO: AO mode (LTAO, SCAO, Gaussian)'
        print '12. res_jitter: Residual telescope jitter [mas]'
        print '13. site_temp: Site/telescope temperature [K]'
        print '14. spec_Nyquist: (True/False) - Set spectral sampling to Nyquist sample spectral resolution element'
        print '15. Spec_samp: Spectral sampling [A/pix] - Only used if spec_Nyquist = False, but enter a value regardless!'
        print '16. noise_force_seed: Force random number seed to take a set value (0=No, 1-10=yes and takes that value)'
        print '17. remove_background: (True/False) - Subtract background spectrum'
        print '18. return_object: (True/False) - Return object cube'
        print '19. return_transmission: (True/False) - Return transmission cube'
        print '20. Turn ADR off: (True/False) - atmospheric differential refraction'
        print ""

    if len(sys.argv) == 21:

        datacube = str(sys.argv[1])
        DIT = int(sys.argv[2])
        NDIT = int(sys.argv[3])
        band = str(sys.argv[4])
        spax = n.array([float(sys.argv[5]), float(sys.argv[6])])
        seeing = float(sys.argv[7])
        zenith_ang = float(sys.argv[8])
        telescope = str(sys.argv[9])
        user_PSF = str(sys.argv[10])
        AO = str(sys.argv[11])
        res_jitter = float(sys.argv[12])
        site_temp = float(sys.argv[13])
        combine_ndits = True
        spec_Nyquist = sys.argv[14]
        if spec_Nyquist == 'True':
            spec_N = True
        elif spec_Nyquist == 'False':
            spec_N = False
        Spec_samp = float(sys.argv[15])
        noise_force_seed = int(float(sys.argv[16]))
        remove_background = sys.argv[17]
        return_obj = sys.argv[18]
        return_tra = sys.argv[19]


        main(datacube, DIT, NDIT, band, spax, seeing, zenith_ang, telescope, user_PSF,
             AO, res_jitter, site_temp, combine_ndits, spec_N, Spec_samp,
             noise_force_seed, remove_background, return_obj, return_tra)
