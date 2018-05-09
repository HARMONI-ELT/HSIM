'''Spectral resolution function

Author: Simon Zieleniewski

Last updated: 18-09-14

'''

import numpy as n
from scipy.ndimage.filters import convolve1d
from modules.Gaussians import Gauss
from modules.frebin import *


def spectral_res(datacube, head, R, band, wavels, spec_nyquist=True, spec_samp=1.):
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
        R: spectral resolving power
        band: string representing wavelength band choice.
              Option of "None" which ignores photometric bands. 
        wavels: wavelength array for datacube
        spec_nyquist: Boolean - if True, 2 pixels per resolution element.
                              - if False, uses spec_samp value

    Outputs:
        scaled_cube: rebinned datacube
        head: updated header file
        new_wavels: new wavelength array for datacube
        delta_lambda: New resolution [um]
        lamb_per_pix: pixel dispersion [um/pixel]
    '''

    print 'Spectral resolution and wavelength range.'
    print 'Chosen resolution: ', R
    print 'Chosen band: ', band

    bands = {'V':(.47,.63),'R':(.63,.79), 'Iz':(.82,1.03),
             'J':(1.08,1.36), 'H':(1.46,1.83), 'K':(1.95,2.45),
             'V+R':(.47,.81), 'Iz+J':(.8,1.36), 'H+K':(1.45,2.45),
             'V-high':(.53,.59), 'R-high':(.61,.68),
             'z':(.82,.91), 'J-high':(1.17,1.29),
             'H-high':(1.545,1.715), 'K-high':(2.09,2.32),
             'None':None}

    bandws = bands[band]

    z, y, x = datacube.shape

    if bandws == None and spec_nyquist==False:
        print 'Ignoring spectral dimension!'
        return datacube, head, wavels, head['SPECRES'], head['CDELT3']
    elif bandws != None and spec_nyquist==True:
        new_res = (bandws[0]+bandws[1])/(2.*R)
        lamb_per_pix = new_res/2.
    elif bandws != None and spec_nyquist==False:
        new_res = (bandws[0]+bandws[1])/(2.*R)
        lamb_per_pix = spec_samp/10000.
    elif bandws == None and spec_nyquist==True:
        #no convolution, just interpolate and read out new wavelengths
        #when sig=0, gauss returns [1.,1.] so convolution does nothing.
#        new_res = head['SPECRES']
        new_res = (wavels[0]+wavels[-1])/(2.*R)
        lamb_per_pix = new_res/2.
    else:
        print 'Ignoring spectral dimension!'
        return datacube, head, wavels, head['SPECRES'], head['CDELT3']        

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

    
    #New wavelength values from chosen photometric band (or ignore)
    if bandws:
        new_wavels = n.arange(bandws[0], bandws[1], lamb_per_pix)
        new_wavels[-1] = bandws[1]
        #Slice datacube down to size of new wavelength array
        start = n.where(wavels < new_wavels[0])[0][-1]
        end = n.where(wavels > new_wavels[-1])[0][0]
        datacube = datacube[start:end+1,:,:]
    elif not bandws:
        new_wavels = n.arange(wavels[0], wavels[-1], lamb_per_pix)
        new_wavels[-1] = wavels[-1]   
        
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
            new_cube[:,:,i] = frebin(datacube[:,:,i], (y,len(new_wavels)),True)

    #print 'New cube sum = ', new_cube.sum()

    #Put datacube flux back into photons/wavelength units
    if head['FUNITS'] == 'erg/s/cm2/A/arcsec2' or head['FUNITS'] == 'J/s/m2/A/arcsec2':
        new_cube /= (lamb_per_pix*10000.)
    else:
        new_cube /= lamb_per_pix
        

    #Update header
    head.update('CRVAL3', new_wavels[0])
    head.update('CDELT3', (new_wavels[1]-new_wavels[0]))
    head.update('NAXIS3', len(new_wavels))
    head.update('SPECRES', new_res)

    print 'Spectral resolution and wavelength range - done!'
    
    #Return new_datacube, new_wavels, updated_header
    return new_cube, head, new_wavels, new_res, lamb_per_pix


