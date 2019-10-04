'''Code to perform low resolution prism mode
for NIRSpec simulation code.

By Simon Zieleniewski

Last updated 01-11-15

'''

import os
import numpy as n
import scipy.interpolate as si
from scipy.ndimage.filters import convolve1d
from modules.Gaussians import Gauss
from modules.misc_utils import path_setup


prism_path = path_setup('../../Sim_data/Prism/')


def nr(lamb, s, Rpower):
    '''Crude function to perform N-R solving method
    on equation: lam - dellam(lam)/2. = s,
    where s = lam_0 + dellam(lam_0)/2,
    and dellam(lam) = lam/R.

    Inputs:

        lamb - wavelength [um]
        s - previous step [um]
        R - Resolving power

    Outputs:

        lam - new wavelength [um]

    '''
    for i in xrange(500):
        lamb = lamb - (2.*(lamb-s) - lamb/float(Rpower))/float(2.-(1/float(Rpower)))
    return lamb



def low_res_mode(cube, head, wavels):
    '''Function to convolve spectral dimension of cube for low resolution, R~100
    mode of NIRSpec, covering wavelength range of 0.6-5.0 um.
    Because low resolution mode has variable resolution, convolution is performed
    in pixel space and cube is returned in pixel space.

    Inputs:

        cube: datacube
        head: datacube header
        wavels: wavelength array [um]

    Outputs:
        scaled_cube: rebinned datacube
        head: updated header file
        new_wavels: new wavelength array for datacube [um]
        delta_lambda: New resolution [um]
        lamb_per_pix: pixel dispersion [um/pixel]
        outR: output Resolving power
    '''

    # - Load pixel wavelength resolution data
    # - Calculate pixel wavelength bins according to input cube resolution and resolution data
    # - cut datacube down in wavelength range if larger than full low res range
    # - if cube shorter than full low res range, process the full wavelength range of cube
    # - bin spectrum up in pixel space
    # - 10x hypersample and convolve with Gaussian LSF in pixel space
    # - resample back to normal sampling in pixel space
    # - Return cube in pixel space, along with header, (irregular) wavelength array, (irregular) delta_lambda array, and
    #(irregular) pixel dispersion array which will be half of each delta_lambda value
    # - output cube, updated header, new wavels, delta_lambda (NEED TO CHECK HOW TO RETURN THIS AS THIS IS USED BY BACKGROUND FUNCTION), lamb_per_pix
        
    #Load data
    dat = n.genfromtxt(os.path.join(prism_path,'Prism_v2.csv'), delimiter=', ')
    dellam = dat[:,0] / dat[:,1]

    #Ignore spectral dimension if input cube is coarser spectral resolution than output
    if head['SPECRES'] > dat[:,2].any():
        print 'Input spectral resolution coarser than output.'
        print 'Ignoring spectral dimension.'       
        return cube, head, wavels, head['SPECRES'], head['CDELT3']
    else:
        convres = n.sqrt(dat[:,2]**2 - head['SPECRES']**2)
        
    outputres = si.interp1d(dat[:,0], dat[:,2])
    rinterp = si.interp1d(dat[:,0], convres)

    #Create new finely sampled wavelength array (0.6-5.0 um) for pixel summation
    newlam = n.arange(0.6, 5.0, dat[0,2]/100.)
    newres = rinterp(newlam)
    newrs = newlam/newres
    rpower = si.interp1d(newlam, newrs)
    dlam = si.interp1d(newlam, newres)

    lams = []
    dlams = []

    lam = newlam[0]
    dl = dlam(lam)
    lams.append(lam)
    dlams.append(dl)
    try:
        for i in xrange(2000):
            lam = nr(lam, lam+dl/2., rpower(lam+dl/2.))
            lams.append(lam)
            dlams.append(dlam(lam))
            dl = dlam(lam)
    except ValueError:
        print 'Pixel wavelength values determined'

    lams = n.array(lams)
    dlams = n.array(dlams)   
    
    #Find relationship between wavelength and pixel no.
    wavetopix = n.column_stack((range(len(lams)), lams))
    wtpinterp = si.interp1d(wavetopix[:,0], wavetopix[:,1])
    #Nyquist sample pixels (2 pixels per lambda)
    nyqpix = n.linspace(wavetopix[0,0], wavetopix[-1,0], 2*len(wavetopix[:,0]))
    wavetonyqpix = n.column_stack((nyqpix, wtpinterp(nyqpix)))
    wtnpinterp = si.interp1d(wavetonyqpix[:,0], wavetonyqpix[:,1])

    #Put datacube flux into photon units
    if head['FUNITS'] == 'erg/s/cm2/A/arcsec2' or head['FUNITS'] == 'J/s/m2/A/arcsec2':
        cube *= (head['CDELT3']*10000.)
    else:
        cube *= head['CDELT3']
    
    #Interpolate datacube onto regular pixel grid (irregular wavelength grid) with 10x hypersampling
    if wavels[0] < wavetonyqpix[0,1] and wavels[-1] > wavetonyqpix[-1,1]:#cube wavelength range larger than low res mode
        print 'Cube wavelength range larger than low res mode'
        start = 0
        end = -1
        incube_interp = si.interp1d(wavels, cube, axis=0)
        hpix = n.linspace(wavetonyqpix[0,0], wavetonyqpix[-1,0], len(wavetonyqpix)*10.)
    elif wavels[0] > wavetonyqpix[0,1] and wavels[-1] < wavetonyqpix[-1,1]:#cube wavelength range inside low res mode
        print 'Cube wavelength range inside low res mode range'
        incube_interp = si.interp1d(wavels, cube, axis=0)
        start = n.where(wavetonyqpix[:,1] > wavels[0])[0][0]
        end = n.where(wavetonyqpix[:,1] > wavels[-1])[0][0] - 1
        hpix = n.linspace(wavetonyqpix[start,0], wavetonyqpix[end,0], len(wavetonyqpix[start:end,0])*10.)
    elif wavels[0] > wavetonyqpix[0,1]:#cube short wavelength longer than low res mode shortest wavelength
        print 'Cube shortest wavelength larger than low res mode'
        incube_interp = si.interp1d(wavels, cube, axis=0)
        start = n.where(wavetonyqpix[:,1] > wavels[0])[0][0]
        end = -1
        hpix = n.linspace(wavetonyqpix[start,0], wavetonyqpix[-1,0], len(wavetonyqpix[start:,0])*10.)
    elif wavels[-1] < wavetonyqpix[-1,1]:#cube longest wavelength shorter than low res mode longest wavelength
        print 'Cube longest wavelength shorter than low res mode'
        incube_interp = si.interp1d(wavels, cube, axis=0)
        start = 0
        end = n.where(wavetonyqpix[:,1] >wavels[-1])[0][0] - 1
        hpix = n.linspace(wavetonyqpix[0,0], wavetonyqpix[end,0], len(wavetonyqpix[:end,0])*10.)
    else:
        raise ValueError('Wavelength error!!!')
   
    hpixlams = wtnpinterp(hpix)   
    hpixcube = incube_interp(hpixlams)

    #Convolve with Gaussian of 2 pix FWHM (hypersampled to 20 hpix FWHM)
    #hgauss = lr_Gauss(2/float(2.*n.sqrt(2.*n.log(2))), 0.1, pix_space=True)
    hgauss = Gauss(2/float(2.*n.sqrt(2.*n.log(2))), 0.1)
    convcube = convolve1d(hpixcube, hgauss[:,1], axis=0)
    
    #Resample back onto normal pixel grid
    hconvcube_interp = si.interp1d(hpix, convcube, axis=0)    
    convcube = hconvcube_interp(wavetonyqpix[start:end,0])

    #final cube
    final_cube = convcube

    #Ensure flux conservation (simplified method)
    print 'Input cube sum = ', n.sum(hpixcube)*n.diff(hpix)[0]
    target = n.sum(hpixcube)*n.diff(hpix)[0]
    current = n.sum(final_cube)*n.diff(wavetonyqpix[start:end,0])[0]
    fac = current/target
    final_cube = n.divide(final_cube, fac)
    print 'Output cube sum = ', n.sum(final_cube)*n.diff(wavetonyqpix[start:end,0])[0]

    #Output wavelength values
    final_lams = wavetonyqpix[start:end,1]

    #Output spectral sampling
    out_samp = outputres(final_lams)/2.
    out_samp.shape = (len(out_samp),1,1)
    
    #Put datacube flux back into photons/wavelength units
    if head['FUNITS'] == 'erg/s/cm2/A/arcsec2' or head['FUNITS'] == 'J/s/m2/A/arcsec2':
        final_cube /= (out_samp*10000.)
    else:
        final_cube /= out_samp

    #Update header
    head.update('CTYPE3', 'NL WAVELENGTH', 'Non-linear wavelength')
    head.update('CRVAL3', final_lams[0])
    head.update('NAXIS3', len(final_lams))
    del head['CDELT3']

    outR = 100
    #Return final cube, updated header, new wavelengths, column_stack of (pixel no. & pixel wavelength values), lambda_per_pixel_array
    return final_cube, head, final_lams, wavetonyqpix[start:end,:], out_samp, outR



def low_res_spec(spec, wavetopix, transmission_spec=False):
    '''Function to take input spectrum (e.g. Skycalc transmission or emission spectrum)
    and return as observed by NIRSpec Prism mode.

    Inputs:

        spec: column stacked spectrum  (wavelength [um], flux)
        wavetopix: Column stacked Nyquist sampling pixels and pixel wavelength values for R=500 mode
        transmission_spec: Boolean. If True, do not multiply then divide by spectral sampling as not required.

    Outputs:

        outspec: Output spectrum (wavelength [um], flux)

    '''

    dat = n.genfromtxt(os.path.join(prism_path,'Prism_v2.csv'), delimiter=', ')
    outputres = si.interp1d(dat[:,0], dat[:,2])

    #Put spectrum flux into photon units (except if transmission spectrum)
    if not transmission_spec:
        specflux = spec[:,1] * (spec[1,0]-spec[0,0])
    elif transmission_spec:
        specflux = spec[:,1]

    
    #Find relationship between wavelength and pixel no.
    wtnpinterp = si.interp1d(wavetopix[:,0], wavetopix[:,1])

    #Interpolate spectrum onto regular pixel grid (irregular wavelength grid) with 10x hypersampling
    if spec[0,0] < wavetopix[0,1] and spec[-1,0] > wavetopix[-1,1]:#cube wavelength range larger than low res mode
        inspec_interp = si.interp1d(spec[:,0], specflux)
#        start = 0
#        end = -1
        hpix = n.linspace(wavetopix[0,0], wavetopix[-1,0], len(wavetopix)*10.)
    elif spec[0,0] > wavetopix[0,1] and spec[-1,0] < wavetopix[-1,1]:#cube wavelength range inside low res mode
        inspec_interp = si.interp1d(spec[:,0], specflux)
        start = n.where(wavetopix[:,1] > spec[0,0])[0][0]
        end = n.where(wavetopix[:,1] > spec[-1,0])[0][0]
        hpix = n.linspace(wavetopix[start,0], wavetopix[end,0], len(wavetopix[start:end,0])*10.)
    elif spec[0,0] > wavetopix[0,1]:#cube short wavelength longer than low res mode shortest wavelength
        inspec_interp = si.interp1d(spec[:,0], specflux)
        start = n.where(wavetopix[:,1] > spec[0,0])[0][0]
#        end = -1
        hpix = n.linspace(wavetopix[start,0], wavetopix[-1,0], len(wavetopix[start:,0])*10.)
    elif spec[-1,0] < wavetopix[-1,1]:#cube longest wavelength shorter than low res mode longest wavelength
        inspec_interp = si.interp1d(spec[:,0], specflux)
#        start = 0
        end = n.where(wavetopix[:,1] > spec[-1,0])[0][0]
        hpix = n.linspace(wavetopix[0,0], wavetopix[end,0], len(wavetopix[:end,0])*10.)
    else:
        raise ValueError('Wavelength error!!!')
    
##    inspec_interp = si.interp1d(spec[:,0], specflux)
##    hpix = n.linspace(wavetopix[0,0], wavetopix[-1,0], len(wavetopix)*10.)
    hpixlams = wtnpinterp(hpix)   
    hpixspec = inspec_interp(hpixlams)

    #Convolve with Gaussian of 2 pix FWHM (hypersampled to 20 hpix FWHM)
    #hgauss = lr_Gauss(2/float(2.*n.sqrt(2.*n.log(2))), 0.1, pix_space=True)
    hgauss = Gauss(2/float(2.*n.sqrt(2.*n.log(2))), 0.1)
    convspec = convolve1d(hpixspec, hgauss[:,1])
    
    #Resample back onto normal pixel grid
    hconvspec_interp = si.interp1d(hpix, convspec)    
    convspec = hconvspec_interp(wavetopix[:,0])

    #final spectrum
    final_spec = n.column_stack((wavetopix[:,1], convspec))
    #Output sampling
    out_samp = outputres(wavetopix[:,1])/2.

    #Ensure flux conservation (simplified method)
    print 'Input spectrum sum = ', n.sum(hpixspec)*n.diff(hpix)[0]
    target = n.sum(hpixspec)*n.diff(hpix)[0]
    current = n.sum(final_spec[:,1])*n.diff(wavetopix[:,0])[0]
    fac = current/target
    final_spec[:,1] = n.divide(final_spec[:,1], fac)
    print 'Output cube sum = ', n.sum(final_spec[:,1])*n.diff(wavetopix[:,0])[0]
    
    #Put spectrum flux back into photons/wavelength units (except transmission spectrum)
    if not transmission_spec:
        final_spec[:,1] /= out_samp
    

    #Return final spectrum
    return final_spec


###Gaussian function to use in low_res_mode
##def lr_Gauss(sigma, delta_x, pix_space=False):
##    '''Function that creates an array of a normalised
##    Gaussian distribution for use in convolutions.
##
##    Inputs:
##        sigma: Gaussian dispersion (related to FWHM by
##               FWHM = 2*sqrt(2ln(2))*sigma
##        delta_x: resolution [x/pixel]
##
##    Outputs:
##        column array of (x, y)
##    '''
##    if sigma != 0.:
##        #Create Gaussian
##        print sigma
##        print sigma*10
##        print float(delta_x)
##        print sigma*10/float(delta_x)
##        num_xs = int(sigma*10/float(delta_x))
##        print num_xs
##        if n.mod(num_xs,2) == 0. and pix_space != True:
##            print 'In if statement'
##            num_xs += 1.
##        print pix_space
##        xs = n.linspace(-5*sigma, 5*sigma, num=num_xs)
##        ys = (1/(n.sqrt(2.*n.pi)*sigma))*n.exp(-0.5*(xs/float(sigma))**2)
##
##        ys /= n.sum(ys)
##
##        return n.column_stack((xs, ys))
##    elif sigma == 0.:
##        return n.array([1.,1.])
