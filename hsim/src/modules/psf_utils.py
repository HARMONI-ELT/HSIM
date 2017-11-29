'''PSF functions

Author: Simon Zieleniewski

Last updated: 28-04-16

'''


import numpy as n
import scipy as s
from frebin import *
from Math_functions import obsc_airy, moffat, lorentz
from Gaussians import Gauss2D


#Constant for converting between radians and milliarcseconds
rad_conv = 1.E-3/206265.


def r_values(spaxels, spans):
    '''Generates an array of size [spaxels,spaxels] containing the distance
       of each element from the centre point of the array.

       =========================== Input ==============================
       spaxels       - List, number of spaxels in [x, y] directions.
                       x-y grid in the format: datacube[Intensity, y, x].
                       (Should be even)
       span          - List, the width of the edges of the image in arc
                       seconds [x, y].
       centre        - BOOL, if True, image is central, if False then
                       image is offset to centre on a pixel at just off
                       the centre of the spaxel array. (currently removed)
       ============================ Out ===============================
       r             - 2D numpy array. Contains the distance of each
                       element from the centre point of the array.'''

    spans = n.array(spans, dtype=n.float32)

    i=(spaxels[0]/2)-1
    j=(spaxels[1]/2)-1
    incrx=spans[0]/(2*spaxels[0])#Half a spaxel width in arcsec
    incry=spans[1]/(2*spaxels[1])#Half a spaxel width in arcsec
    x=n.linspace(-(spans[0]/2.),i*2*incrx,spaxels[0])
    y=n.linspace(-(spans[1]/2.),i*2*incry,spaxels[1])

    cart=n.meshgrid(x,y)
    r=n.sqrt(cart[0]**2+cart[1]**2)    #Pythagoras
    r /= 206265.                        #Values in radians
    #Quick and dirty fix to ensure central value is tiny for Airy function to handle.
    r = n.where(r < 1.E-15, 1.E-30, r)

    return r



def dist_arr(size, spax):
    '''Function that generates array with distances from centre

    Inputs:

        - size: length of array in pixels
        - spax: spaxel size [mas/spaxel]

    Outputs:

        - distance_array: shape (size,size)

    '''

    arra = n.zeros((size,size),dtype=float)
    y, x = arra.shape
    x_inds, y_inds = n.ogrid[:y, :x]
    mid_x, mid_y = (s.array(arra.shape[:2]) - 1) / float(2)
    n_arra = ((y_inds - mid_y) ** 2 + (x_inds - mid_x) ** 2) ** 0.5
    n_arra *= (spax/206265.)
    #Quick and dirty fix to ensure central value is tiny for Airy function to handle.
    r = n.where(n_arra < 1.E-15, 1.E-30, n_arra)
    return r



def residual_jitter(array, array_spax, res_jitt):
    '''Function that takes wavelength channel array and convolves with
    Gaussian of FWHM given by res_jitt.

    Inputs:

        array: 2D array corresponding to wavelength channel
        array_spax: Spatial scale of array [mas/pix]
        res_jitt: Residual jitter sigma in units of [mas]
                  FWHM = 2.*n.sqrt(2.*n.log(2.))*sigma

    Outputs:

        jitt_array: New array of jittered PSF.
    '''

    #If res_jitt == 0, ignore function and return array.
    if res_jitt == 0.:
        return array
    else:
        fwhm_jitt = 2.*n.sqrt(2.*n.log(2.)) * res_jitt
        pix_fwhm = abs(fwhm_jitt)/(array_spax)
        g_array = Gauss2D(array.shape[0], pix_fwhm)

        #Convolve using FFTs
        conv_array = n.fft.ifft2(n.fft.fft2(array)*n.fft.fft2(g_array))
        jitt_array = n.fft.fftshift(conv_array.real)

        return jitt_array



def create_Gausspsf_channel(wave, seeing, aperture, res_jitter, array_dim, Nyquist, psfspax, output_spax):
    '''Function that creates a 2D spatial Gaussian PSF.
    The FWHM of the Gaussian is depends on wavelength in the following way:

    FWHM(radians) = 0.98 Lambda/r_0
    r_0 = A Lambda^(6/5.)
    FWHM(arcsec) = (0.98/A)/Lambda^(1/5.)

    Inputs:

        wave: wavelength [microns]
        seeing: FWHM value of the seeing at 500nm (V-band) in arcseconds.
        aperture: List containing [diameter of telescope, obscuration ratio].
                  Diameter only used when Nyquist=True.
        array_dim: Spatial size of PSF arrays.
        Nyquist: Boolean - True returns Nyquist sampling for each wavelength channel.
                         - False uses value of spaxel.
        psfspax: Initial spatial sampling value [mas].
        output_spax: Output spatial sampling [mas].


    Output:

        psf: PSF array for given wavelength.

    '''

##    lam = 500.E-9
##    #Compute r_0 from given seeing value
##    r_0 = 0.98*lam/float(seeing/206265.)
##
##    #Use r_0 to find A
##    A = r_0/lam**(6/5.)
##
##    #Use A to find constant for FWHM propto Lambda**(-1/5) equation
##    B = 0.98*206265./A

    r0 = 0.976*500.E-9/float(seeing)*(180./n.pi*3600.)*(wave/0.5)**1.2 #in m
    fwhm = seeing*(wave/0.5)**(-0.2)*n.sqrt(1.+(1./(1.+300.*aperture[0]/23.)-1)*2.183*(r0/23.)**(0.356)) #in arcsec

    if Nyquist:
        #fwhm = B/float((wave*(1.E-6))**(1/5.)) #in arcsec
        diff = wave*(1.E-6)/(2.*float(aperture[0]))
        diff *= 206265 #in arcsec/spaxel
        channel = Gauss2D(array_dim, fwhm/diff)
        #Apply residual jitter
        psf = residual_jitter(channel, diff, res_jitter)

    elif not Nyquist:
        spax = psfspax / 1000. #in arcsec/spaxel
        #fwhm = B/float((wave*(1.E-6))**(1/5.)) #in arcsec
        channel = Gauss2D(array_dim, fwhm/spax)
        #Apply residual jitter
        psf = residual_jitter(channel, spax, res_jitter)

    #Normalise PSF
    psf /= psf.sum()

    #Frebin PSF up to same sampling as datacube channels

    #total field of view in mas
    x_field = array_dim*psfspax
    y_field = array_dim*psfspax
    x_newsize = x_field/float(output_spax[0])
    y_newsize = y_field/float(output_spax[1])

    psf = frebin(psf, (x_newsize, y_newsize), total=True)

    return psf



def create_psf_channel(params, iteration, psfspax, output_spax, array_dim, aperture, res_jitter):
    '''Function that creates an AO PSF image given the
    parameter values and sampling scale.

    Inputs:

        psfparams: Dictionary containing parameter arrays
        iteration: Iteration value through wavelength array
        psfspax: PSF sampling scale [mas/spaxel]
        output_spax: Output sampling scale tuple (x, y) [mas/spaxel]
        array_dim: Size of PSF array
        aperture: List containing telescope diameter [m] and obsc. ratio
        res_jitter: value of residual telescope jitter [mas]

    Outputs:

        psf = 2D PSF array

    '''

    #array = r_values([array_dim, array_dim], [(psfspax/1000.)*array_dim,(psfspax/1000.)*array_dim])
    array = dist_arr(array_dim, psfspax/1000.)

    airy = obsc_airy(array, params['oh'][iteration], params['ow'][iteration]*rad_conv, aperture[1])
    moff1 = moffat(array, params['mh'][iteration], params['mw'][iteration]*rad_conv, params['mq'][iteration])
    lor = lorentz(array, params['lh'][iteration], params['lp'][iteration]*rad_conv, params['lw'][iteration]*rad_conv)
    moff2 = moffat(array, params['m2h'][iteration], params['m2w'][iteration]*rad_conv, params['m2q'][iteration])
    channel = (airy + moff1 + lor + moff2)
    #Apply residual jitter
    channel /= channel.sum() #normalise
    psf = residual_jitter(channel, psfspax, res_jitter)

    #Frebin PSF up to same sampling as datacube channels
    psf /= psf.sum() #normalise

    #total field of view in mas
    x_field = array_dim*psfspax
    y_field = array_dim*psfspax
    x_newsize = x_field/float(output_spax[0])
    y_newsize = y_field/float(output_spax[1])

    psf = frebin(psf, (x_newsize, y_newsize), total=True)

    return psf



def create_instpsf(array_dim, psfspax, output_spax, psfoutspax):
    '''Function that creates an instrument PSF array given the
    parameter values and sampling scale.

    Inputs:

        array_dim: Size of PSF array
        psfspax: PSF sampling scale [mas/spaxel]
        output_spax: Output sampling scale tuple (x, y) [mas/spaxel]
        psfoutspax: PSF convolution sampling scale to match datacube: tuple (x, y) [mas/spaxel]

    Outputs:

        instpsf = 2D PSF array

    '''

    #FWHM of Instrument PSF depending on output spaxel scale
    #Factors taken into account:
    #design image quality, manufacturing and assembly tolerances, vibration, flexure, diff refraction,
    #Also interpolating the data back onto a regular grid.
    #The interpolation adds another contribution of 0.8-1.8 PIXEL FWHM.
    #So a reasonable assumption: 1.1 pix FWHM.

    #FWHMs given in mas according to output spaxel scale
    instpsf_data = {(30.,60.): 65.,
                    (20.,20.): 34.,
                    (10.,10.):17.,
                    (4.,4.):6.}
    #print 'SPAX = ', output_spax
    try:
        FWHM = instpsf_data[output_spax]
        #print 'FWHM = ', FWHM
        instpsf = Gauss2D(array_dim, FWHM/float(psfspax))
    except:
        #print 'Non-HARMONI spaxel scale - no instrument PSF'
        instpsf = n.zeros([array_dim, array_dim], dtype=float)
        instpsf[array_dim/2,array_dim/2] = 1.0

    #total field of view in mas
    x_field = array_dim*psfspax
    y_field = array_dim*psfspax
    x_newsize = x_field/float(psfoutspax[0])
    y_newsize = y_field/float(psfoutspax[1])

    instpsf = frebin(instpsf, (x_newsize, y_newsize), total=True)

    return instpsf
