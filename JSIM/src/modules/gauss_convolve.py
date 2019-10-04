'''Gaussian convolve

Author: Simon Zieleniewski

Last updated: 29-09-14

'''

import numpy as n
import scipy as s
from Gaussians import Gauss


def gauss_convolve(SED, sigma, lambda_space='Linear'):
    '''Function that convolves an array with a normalised Gaussian
    of width 2*sqrt(2*ln(2))*sigma. Option of either linear or log spaced
    convolution.

    Inputs:
        SED: column array of wavelength and corresponding value.
        sigma: sigma value of Gaussian function.
        lambda_space: Two options - 'Linear' or 'Log'. Sets
                      convolution space.

    Outputs:
        new_SED: new convolved column array. Always returned
                 in linear space.
        '''

    #Interpolate SED
    interspec = s.interpolate.interp1d(SED[:,0], SED[:,1], kind='linear')

    if lambda_space=='Linear':
        #Get linear spaced wavelength array
        lin_lambs = n.linspace(SED[0,0], SED[-1,0], len(SED[:,0]))
        #Create normalised Gaussian
        gauss_array = Gauss(sigma, (lin_lambs[1]-lin_lambs[0]))
        #Flux values for linear spacing
        lin_ys = interspec(lin_lambs)
        #Convolve linearly spaced flux values with Gaussian
        conv_ys = n.convolve(lin_ys, gauss_array[:,1], mode='same')
        
        return n.column_stack((lin_lambs, conv_ys))
    
    elif lambda_space=='Log':
        #Get log spaced wavelength array
        log_lambs = n.linspace(n.log10(SED[0,0]), n.log10(SED[-1,0]), len(SED[:,0]))
        lambs = 10.**log_lambs
        lambs = n.round(10**log_lambs, decimals=3)
        #Flux values for log spacing
        log_ys = interspec(lambs)
        #Create normalised Gaussian
        gauss_array = Gauss(sigma, (log_lambs[1]-log_lambs[0]))
        #Convolve
        conv_ys = n.convolve(log_ys, gauss_array[:,1], mode='same')
        #Return to linearly spaced array
        lin_lambs = n.linspace(lambs[0], lambs[-1], len(SED[:,0]))
        interspec2 = s.interpolate.interp1d(lambs, conv_ys, kind='linear')
        final_ys = interspec2(lin_lambs)
        
        return n.column_stack((lin_lambs, final_ys))

