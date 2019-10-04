'''Spectral resolution function

Author: Simon Zieleniewski

Last updated: 01-12-15

'''

import numpy as n
from scipy.ndimage.filters import convolve1d
from modules.Gaussians import Gauss
from modules.frebin import *


def spectral_res(datacube, head, grating, wavels, gratings, spec_nyquist=True, spec_samp=1., ignoreLSF=False):
    '''Function that takes input datacube and rebins it to the
    chosen spectral resolution. It convolves datacube along wavelength
    axis, interpolates all spaxels along wavelength axis then extracts
    pixel vales for each datacube wavelength channel. Function combines both
    spectral resolution choice and datacube chopping in wavelength axis.
    Option of bypassing spectral processing by setting band=='None' and
    spec_nyquist=True - this ignores wavelength dimension.

    Inputs:
        datacube: object datacube
        head: header file for object datacube
        grating: string representing chosen grating. Grating choice sets R and wavelength range
        wavels: wavelength array for datacube
        gratings: dictionary of grating parameters (lambda_min, lambda_max, R)
        spec_nyquist: Boolean - if True, 2 pixels per resolution element.
                              - if False, uses spec_samp value
        ignoreLSF: Boolean - ignore LSF convolution in spectral dimension

    Outputs:
        scaled_cube: rebinned datacube
        head: updated header file
        new_wavels: new wavelength array for datacube
        delta_lambda: New resolution [um]
        lamb_per_pix: pixel dispersion [um/pixel]
        resol: output R
    '''

    print 'Spectral resolution and wavelength range.'
    print 'Chosen grating: ', grating

    bandws = gratings[grating]

    z, y, x = datacube.shape

    if ignoreLSF==True:
        print 'Ignoring spectral dimension!'
        return datacube, head, wavels, head['SPECRES'], head['CDELT3'], 0
    if bandws == None and spec_nyquist==False:
        print 'Ignoring spectral dimension!'
        return datacube, head, wavels, head['SPECRES'], head['CDELT3'], 0
    elif bandws != None and spec_nyquist==True:
        #new_res = (bandws[0]+bandws[1])/(2.*bandws[2])
        new_res = bandws[2]
        lamb_per_pix = new_res/2.
    elif bandws != None and spec_nyquist==False:
        #new_res = (bandws[0]+bandws[1])/(2.*bandws[2])
        new_res = bandws[2]
        lamb_per_pix = spec_samp/10000.
    elif bandws == None and spec_nyquist==True:
        print 'Ignoring spectral dimension!'
        return datacube, head, wavels, head['SPECRES'], head['CDELT3'], 0
    else:
        print 'Ignoring spectral dimension!'
        return datacube, head, wavels, head['SPECRES'], head['CDELT3'], 0

    print 'Current resolution = %.5f microns' % head['SPECRES']
    print 'Current sampling = %.5f microns' % head['CDELT3']
    print 'New resolution = %.5f microns' % new_res
    print 'New sampling = %.5f microns' % lamb_per_pix

    if head['SPECRES'] >= new_res:
        print 'WARNING: Input spectral resolution is coarser (or equal) than chosen output!'
        #Chose whether to convolve with LSF or simply resample.
        condition = raw_input('Do you want to convolve with Gaussian LSF of FWHM'
                                  ' given by band and chosen resolving power [y], or not [any other key]?: ')

        if condition=='y':
            print 'Convolving with LSF'
            #Generate Gaussian LSF of FWHM = (bandws[0]+bandws[1])/(2.*R)
            sig = new_res/(2.*n.sqrt(2.*n.log(2.)))
            gauss_array = Gauss(sig, head['CDELT3'])
            #Convolve datacube array with Gaussian along wavelength axis
            datacube = convolve1d(datacube, gauss_array[:,1], axis=0)

    else:
        print 'Input spectral resolution smaller than chosen output.'
        print 'Convolving with corresponding LSF.'
        #Generate Gaussian LSF of FWHM = sqrt(new_resolution**2-head['SPECRES']**2)
        sig = n.sqrt(new_res**2-head['SPECRES']**2)/(2.*n.sqrt(2.*n.log(2.)))
        gauss_array = Gauss(sig, head['CDELT3'])

        #Convolve datacube array with Gaussian along wavelength axis
        datacube = convolve1d(datacube, gauss_array[:,1], axis=0)


##    #New wavelength values from chosen photometric band (or ignore)
##    if bandws:
##        new_wavels = n.arange(bandws[0], bandws[1], lamb_per_pix)
##        new_wavels[-1] = bandws[1]
##        #Slice datacube down to size of new wavelength array
##        start = n.where(wavels < new_wavels[0])[0][-1]
##        end = n.where(wavels > new_wavels[-1])[0][0]
##        datacube = datacube[start:end+1,:,:]
##    elif not bandws:
##        new_wavels = n.arange(wavels[0], wavels[-1], lamb_per_pix)
##        new_wavels[-1] = wavels[-1]

    #####
    #Interpolate datacube onto regular pixel grid (irregular wavelength grid) with 10x hypersampling
    if wavels[0] <= bandws[0] and wavels[-1] >= bandws[1]:#cube wavelength range larger than grating choice
        print 'Cube wavelength range larger than grating: chopping down to grating range'

        new_wavels = n.arange(bandws[0], bandws[1]+lamb_per_pix/2., lamb_per_pix)
        new_wavels[-1] = bandws[1]
        start = n.where(wavels <= new_wavels[0])[0][-1]
        end = n.where(wavels >= new_wavels[-1])[0][0]
        #datacube = datacube[start:end+1,:,:]

    elif wavels[0] > bandws[0] and wavels[-1] < bandws[1]:#cube wavelength range inside grating choice
        print 'Cube wavelength range inside grating range'

        new_wavels = n.arange(wavels[0], wavels[-1]+lamb_per_pix/2., lamb_per_pix)
        new_wavels[-1] = wavels[-1]
        start = n.where(bandws[0] < wavels)[0][0]
        end = n.where(bandws[1] > wavels)[0][-1]
        #datacube = datacube[start:end+1,:,:]

    elif wavels[0] > bandws[0]:#cube short wavelength longer than grating shortest wavelength
        print 'Cube shortest wavelength larger than grating'

        new_wavels = n.arange(wavels[0], bandws[1]+lamb_per_pix/2., lamb_per_pix)
        new_wavels[-1] = bandws[1]
        start = n.where(bandws[0] < wavels)[0][0]
        end = n.where(bandws[1] > wavels)[0][-1]
        #datacube = datacube[start:end+1,:,:]

    elif wavels[-1] < bandws[1]:#cube longest wavelength shorter than grating longest wavelength
        print 'Cube longest wavelength shorter than grating'

        new_wavels = n.arange(bandws[0], wavels[-1]+lamb_per_pix/2., lamb_per_pix)
        new_wavels[-1] = wavels[-1]
        start = n.where(bandws[0] < wavels)[0][0]
        end = n.where(bandws[1] > wavels)[0][-1]
        #datacube = datacube[start:end+1,:,:]

    else:
        raise ValueError('Specres wavelength error!!!')

    #####

    new_cube = n.zeros((len(new_wavels), y, x), dtype=float)

    ###New 30-05-14
    #Put datacube flux into photon units
    if head['FUNITS'] == 'erg/s/cm2/A/arcsec2' or head['FUNITS'] == 'J/s/m2/A/arcsec2':
        datacube *= (head['CDELT3']*10000.)
    else:
        datacube *= head['CDELT3']

    #print 'Datacube sum = ', datacube.sum()

    #Bin up in spectral dimension to new pixel sampling
    #Iterate over x-axis and use frebin function on 2D y-z arrays
    for i in xrange(x):
            new_cube[:,:,i] = frebin(datacube[start:end+1,:,i], (y,len(new_wavels)), True)

    #print 'New cube sum = ', new_cube.sum()

    #Put datacube flux back into photons/wavelength units
    if head['FUNITS'] == 'erg/s/cm2/A/arcsec2' or head['FUNITS'] == 'J/s/m2/A/arcsec2':
        new_cube /= (lamb_per_pix*10000.)
    else:
        new_cube /= lamb_per_pix


    #Update header
    head['CRVAL3'] = new_wavels[0]
    head['CDELT3'] = (new_wavels[1]-new_wavels[0])
    head['NAXIS3'] = len(new_wavels)
    head['SPECRES'] = new_res

    resol = n.round(new_wavels[len(new_wavels)/2]/bandws[2], 0)

    print 'Spectral resolution and wavelength range - done!'

    #Return new_datacube, new_wavels, updated_header
    return new_cube, head, new_wavels, new_res, lamb_per_pix, resol
