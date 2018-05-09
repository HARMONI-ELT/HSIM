'''Module that contains the loop over wavelength channels for the
simulation pipeline.

Written by Simon Zieleniewski

Last updated 19-01-17

'''

import numpy as n
from scipy.interpolate import interp1d
from TheSimulator_ADR import add_ADR
from modules.frebin import *
from modules.Gaussians import Gauss2D
from modules.misc_utils import update_progress
from modules.psf_utils import create_Gausspsf_channel, create_psf_channel, create_instpsf
from modules.psf_convolution import psf_convolve
import multiprocessing as mp


def wavelength_loop(cube, head, wavels, out_cube, AO, psfvars, adrvals, newsize, outspax, adr_switch='ON'):
    '''Function to take input datacube and process it iteratively through
    each wavelength channel as follows:
    - Generate PSF array for given channel
    - Add effect of ADR to channel
    - Convolve cube channel with PSF
    - Frebin up to chosen output spaxel scale

    Inputs:

        cube: Datacube
        head: Datacube header
        wavels: Wavelength array
        out_cube: Empty output cube
        AO: AO mode [LTAO, SCAO, Gaussian]
        psdvars: list containing [psfparams, psfspax, psfsize, [D,eps],
                                  res_jitter, seeing, user_psf, user_psflams]
        adrvals: array of ADR values
        newsize: tuple containing (x_newsize, y_newsize) array sizes
        outspax: tuple containing (x, y) output spaxel scales
        adr_switch: ON or OFF. Turns ADR effect on or off.

    Output:

        cube: Processed cube
        head: Updated header
        inspax: Input spaxel scale (mas, mas)
        outspax: Output spaxel scale (mas, mas)

    '''

    for l in xrange(len(wavels)):
        #Print percentage bar
        update_progress(n.round(l/float(len(wavels)),2))

        #Create PSF channel
        #If user PSF and 2D image
        if psfvars[-2] != 'None' and psfvars[-1] == 'None':
            psf = psfvars[-2]
        #If User PSF and 3D cube
        elif psfvars[-2] != 'None' and psfvars[-1] != 'None':
            #Interpolate PSF cube
            interp = interp1d(psfvars[-1], psfvars[-2], axis=0)
            psf = interp(wavels[l])

        elif AO == 'LTAO' or AO == 'SCAO' or AO == 'GLAO':
            psf = create_psf_channel(psfvars[0], l, psfvars[1], (head['CDELT1'],head['CDELT2']),
                                     psfvars[2], psfvars[3], psfvars[4])

        elif AO == 'Gaussian':
            psf = create_Gausspsf_channel(wavels[l], psfvars[5], psfvars[3], psfvars[4], psfvars[2], Nyquist=False,
			   psfspax=psfvars[1], output_spax=(head['CDELT1'],head['CDELT2']))

        else:
            print 'AO = ', AO
            raise ValueError('AO choice error!')

        #Create instrument PSF array
        instpsf = create_instpsf(psfvars[2], psfvars[1], outspax, (head['CDELT1'],head['CDELT2']))

        #Add ADR effect to channel
        if adr_switch == 'ON':
            cube[l,:,:] = add_ADR(cube[l,:,:], head, adrvals[l], 'spline3')

        #Convolve cube channel with PSF channel
        channel = psf_convolve(cube[l,:,:], psf)

        #Convolve cube channel with instrument PSF
        channel = psf_convolve(channel, instpsf)

        #Frebin datacube up to output spaxel scale
        newsize = (int(newsize[0]), int(newsize[1]))
        channel *= (head['CDELT1']*head['CDELT2']*1.E-6)
        channel = frebin(channel, (newsize[0],newsize[1]), total=True)
        channel /= (outspax[0]*outspax[1]*1.E-6)

        #Correct ADR effect
        if adr_switch == 'ON':
            adrhead = head.copy()
            adrhead['CDELT1'] = outspax[0]; adrhead['CDELT2'] = outspax[1]
            channel = add_ADR(channel, adrhead, -1.*adrvals[l], 'spline3')

        #Add channel to output cube
        out_cube[l,:,:] = channel

    return out_cube, head



def pp_wavelength_channel(chann, head, wave, l, AO, psfvars, adrval, newsize, outspax, adr_switch):
    '''Function to take input datacube and process it iteratively through
    each wavelength channel as follows:
    - Generate PSF array for given channel
    - Add effect of ADR to channel
    - Convolve cube channel with PSF
    - Frebin up to chosen output spaxel scale

    Inputs:

        chann: cube channel
        head: Datacube header
        wave: wavelength [um]
        l: iteration no.
        out_cube: Empty output cube
        AO: AO mode [LTAO, SCAO, Gaussian]
        psfvars: list containing [psfparams, psfspax, psfsize,
                                  [D,eps], res_jitter, seeing, user_psf,
                                  user_psflams]
        adrval: ADR value
        newsize: tuple containing (x_newsize, y_newsize) array sizes
        outspax: tuple containing (x, y) output spaxel scales
        adr_switch: On or OFF. Turns ADR effect on or off

    Output:

        cube: Processed cube
        head: Updated header
        inspax: Input spaxel scale (mas, mas)
        outspax: Output spaxel scale (mas, mas)

    '''

    #Create PSF channel
    #If user PSF and 2D image
    upsf = psfvars[-2]
    upsflams = psfvars[-1]
    if upsf != 'None' and upsflams == 'None':
        psf = upsf
    #If User PSF and 3D cube
    elif upsf != 'None' and upsflams != 'None':
        #Interpolate PSF cube
        interp = interp1d(upsflams, upsf, axis=0)
        psf = interp(wave)

    elif AO == 'LTAO' or AO == 'SCAO' or AO == 'GLAO':
        psf = create_psf_channel(psfvars[0], l, psfvars[1], (head['CDELT1'],head['CDELT2']),
                                 psfvars[2], psfvars[3], psfvars[4])

    elif AO == 'Gaussian':
        psf = create_Gausspsf_channel(wave, psfvars[5], psfvars[3], psfvars[4], psfvars[2], False,
                       psfvars[1], (head['CDELT1'],head['CDELT2']))

    else:
        print 'AO = ', AO
        raise ValueError('AO choice or user_PSF error!')

    #Create instrument PSF array
    instpsf = create_instpsf(psfvars[2], psfvars[1], outspax, (head['CDELT1'],head['CDELT2']))

    #Add ADR effect to channel
    if adr_switch == 'ON':
        chann = add_ADR(chann, head, adrval, 'spline3')

    #Convolve cube channel with PSF channel
    chann = psf_convolve(chann, psf)

    #Convolve cube channel with instrument PSF
    chann = psf_convolve(chann, instpsf)

    #Frebin datacube up to output spaxel scale
    newsize = (int(newsize[0]), int(newsize[1]))
    chann *= (head['CDELT1']*head['CDELT2']*1.E-6)
    chann = frebin(chann, (newsize[0],newsize[1]), total=True)
    chann /= (outspax[0]*outspax[1]*1.E-6)

    #"Correct" ADR effect
    if adr_switch == 'ON':
        adrhead = head.copy()
        adrhead['CDELT1'] = outspax[0]; adrhead['CDELT2'] = outspax[1]
        chann = add_ADR(chann, adrhead, -1.*adrval, 'spline3')

    return chann, l



def pp_wavelength_loop(cube, head, wavels, out_cube, AO, psfvars, adrvals, newsize, outspax, adr_switch='ON', usecpus=mp.cpu_count()-1):
    '''Function to take input datacube and process it iteratively through
    each wavelength channel as follows:
    - Generate PSF array for given channel
    - Add effect of ADR to channel
    - Convolve cube channel with PSF
    - Frebin up to chosen output spaxel scale

    Inputs:

        cube: Datacube
        head: Datacube header
        wavels: Wavelength array
        out_cube: Empty output cube
        AO: AO mode [LTAO, SCAO, Gaussian]
        psdvars: list containing [psfparams, psfspax, psfsize,
                                  [D,eps], res_jitter, seeing]
        adrvals: array of ADR values
        newsize: tuple containing (x_newsize, y_newsize) array sizes
        outspax: tuple containing (x, y) output spaxel scales
        adr_switch: On or OFF. Turns ADR effect on or off
        usecpus: no. of CPUs to use

    Output:

        cube: Processed cube
        head: Updated header
        inspax: Input spaxel scale (mas, mas)
        outspax: Output spaxel scale (mas, mas)

    '''
    import pprocess as pp

    print 'Using ', usecpus, ' CPUs'
    queue = pp.Queue(limit=usecpus)
    waveloop = queue.manage(pp.MakeParallel(pp_wavelength_channel))

    print "Initialising..."
    for lam in xrange(len(wavels)):
        waveloop(cube[lam,:,:], head, wavels[lam], lam, AO, psfvars, adrvals[lam], newsize, outspax, adr_switch)

    print "Processing..."
    for chan, it in queue:
        update_progress(n.round(it/float(len(wavels)),2))
        out_cube[it,:,:] = chan

    return out_cube, head




if __name__=="__main__":

    import sys
    import time
    import astropy.io.fits as p
    import numpy as n
    from modules.fits_utils import *#fits_header_check, wavelength_array
    from TheSimulator_ADR import *#calc_ref, calc_adr
    from TheSimulator_PSFs import generate_psfcube#generate_psfcube, Gausspsf, r_values

    cube, head = p.getdata(sys.argv[1], header=True)
    fits_header_check(head)
    wavels, head = wavelength_array(head)
    print 'Input cube shape = ', cube.shape

    AO = 'LTAO'
    seeing = 0.67
    zenith_ang = 30.
    D = 37.
    eps = 0.3
    res_jitter = 0.0
    psfspax = 1.
    psfparams, psfsize = generate_psfcube(wavels, AO, seeing, zenith_ang, [D, eps], res_jitter,
                                   head, Nyquist=False, samp=psfspax)
    psfvars = [psfparams, psfspax, psfsize, [D, eps], res_jitter, seeing]

    ppwv = 17.0
    press = 760.0
    site_temp = 280.5
    airmass = 1./n.cos(n.radians(zenith_ang))
    refr = calc_ref(wavels, site_temp, ppwv, press)
    adr = calc_adr(wavels, refr, airmass)

    spax = n.array([10.,10.])
    x_field = head['NAXIS1']*head['CDELT1']
    y_field = head['NAXIS2']*head['CDELT2']
    x_newsize = x_field/float(spax[0])
    y_newsize = y_field/float(spax[1])



    out_cube = n.zeros((len(wavels),int(head['CDELT2']*head['NAXIS2']/spax[1]),
                        int(head['CDELT1']*head['NAXIS1']/spax[0])), dtype=n.float64)
    print 'Output cube shape = ', out_cube.shape


    #Linear operation
    print "Linear operation..."
    start_time = time.time()
    out_cube, head = wavelength_loop(cube, head, wavels, out_cube, AO, psfvars, adr, (x_newsize, y_newsize), spax)

    head.update('CDELT1', spax[0], "mas")
    head.update('CDELT2', spax[1], "mas")
    head.update('NAXIS1', int(x_newsize))
    head.update('NAXIS2', int(y_newsize))

    print "Time elapsed: ", time.time() - start_time, "s"

    p.writeto('./Output_cubes/LINEAR_TEST_CUBE.FITS', out_cube, head, clobber=True)


    #parallel operation
    cube, head = p.getdata(sys.argv[1], header=True)
    fits_header_check(head)
    wavels, head = wavelength_array(head)
    out_cube = n.zeros((len(wavels),int(head['CDELT2']*head['NAXIS2']/spax[1]),
                        int(head['CDELT1']*head['NAXIS1']/spax[0])), dtype=n.float64)
    print "Parallel operation..."
    start_time = time.time()
    out_cube, head = pp_wavelength_loop(cube, head, wavels, out_cube, AO, psfvars, adr, (x_newsize, y_newsize), spax)

    head.update('CDELT1', spax[0], "mas")
    head.update('CDELT2', spax[1], "mas")
    head.update('NAXIS1', int(x_newsize))
    head.update('NAXIS2', int(y_newsize))

    print "Time elapsed: ", time.time() - start_time, "s"

    p.writeto('./Output_cubes/PP_TEST_CUBE.FITS', out_cube, head, clobber=True)
