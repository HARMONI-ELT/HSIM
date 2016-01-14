'''Module to calculate and apply atmospheric differential refraction
to datacubes.

Written by Sarah Kendrew & Simon Zieleniewski

Last updated 31-07-15

'''

import numpy as np
import scipy.ndimage.interpolation as si
import os
import pdb
import math


#====================================================
def calc_ref(l, T, ppwv=17.0, P=760.0):
    
    '''This function calculates the refractive index of air vs. wavelength.
    INPUTS:
    l:      wavelength axis (microns)
    T:      temperature (Kelvin)
    ppwv:    partial pressure of water (millibar)
    P:      atmospheric pressure (millibar)
    
    OUTPUT:
    refractive index of air on the provided wavelength scale'''
    
    Ps = 1013.25
    Ts = 288.15
    
    nout = 1 + ( 64.328 + (29498.1)/(146.0 - l**(-2)) + (255.4)/(41 - l**(-2)) )\
           * ((P*Ts)/(Ps*T)) * 1e-6 - 43.49 * (1.0 - (7.956e-3)/l**2)*(ppwv/Ps) * 1e-6
    
    return nout
#====================================================

def calc_adr(wave, n, am, optlam):
    
    ''' calculates the differential refration in MILLI-ARCSEC
    across the wavelength range provided given the refractive
    index and airmass.
    
    INPUTS:
    wave:       wavelength axis, in microns
    n:          refractive index, unitless (must have same length as wave)
    am:         airmass, unitless (= sec z)
    optlam:     guiding wavelength, in microns
    '''


##    # identify the centre of the band (find the median value)
##    wavec = np.median(wave)
##    wavearg = np.where(wave > wavec)[0][0]
##    nc = n[wavearg]
##    print 'Guiding wavelength: {0}' .format(wavec)

    wavearg = np.where(wave > optlam)[0][0]
    nc = n[wavearg]
    print 'Guiding wavelegth: %.3f' % optlam

    z = np.arccos(1./am)
    
    diffref = 206265. * 1.e3 * ((n**2-1.0)/(2.0*n**2) - (nc**2-1.0)/(2.*nc**2)) * np.tan(z)
    
    return diffref



def add_ADR(array, head, adr_val, interp_mode='spline3'):
    '''Function that takes datacube and applies the effects of
    atmospheric differential refraction by shifting wavelength
    channels along Declination axis.

    Inputs:

        array: Datacube wavelength channel (DEC, RA)
        head: FITS header for datacube
        adr_val: value of ADR shift value for channel
                 relative to the central channel [milliarcsec]
        interp_mode: Choice of 'linear', 'spline3', 'FFT'

    Outputs:

        adr_array: channel with ADR shift

    '''

    
    #check units of both axes
    if head['CUNIT1'] == 'mas' and head['CUNIT2'] == 'mas':
        u_conv = 1.
    elif head['CUNIT1'] == 'arcsec' and head['CUNIT2'] == 'arcsec':
        u_conv = 1.E3
    else:
        raise ValueError('Units of spatial axes are different!')
        

    #Find longest spatial axis to disperse along
    #Displacement in pixel units
    if head['NAXIS1'] >= head['NAXIS2']:
        xdisp = adr_val/float(head['CDELT1'])
        ydisp = 0.0
    else:
        xdisp = 0.0
        ydisp = adr_val/float(head['CDELT2'])

    #Apply shift through interpolation
    if interp_mode == 'linear':
         omode = 1
    elif interp_mode == 'spline3':
        omode = 3
    elif interp_mode == 'FFT':
        #Do something
        pass
    else:
        raise ValueError('Please choose from "linear", "spline3" or "FFT"!')

    adr_channel = si.shift(array*u_conv, [ydisp, xdisp], order=omode, mode='constant', cval=0.0)

    return adr_channel



def optimalguide(wave0, wave1, temp) :
    '''
    Calculates the optimal guiding wavelength for to given end wavelengths
    '''

    # Only need to work with refractive index, as that's the only wavelength dependant term

    nLow = calc_ref(wave0, temp)
    nHigh = calc_ref(wave1, temp)
    nMid = nLow + (nHigh-nLow)/2.0

    # iterate to find the optimal wavelength

    diff=0.0
    dw = wave1-wave0
    waveg = wave0 + (wave1-wave0)/2.0

    while (True) : # 1% tolerance on finding where nMid is??
            waveg=min(max(wave0,waveg + min(1.0,dw)*diff*100000.0),wave1)
            nwave = calc_ref(waveg, temp)
            diff=nwave-nMid
    #		print waveg,nwave,nMid,diff,(nHigh-nLow)*0.001
            if(math.fabs(diff)<(math.fabs(nHigh-nLow))*0.001) : break

    #print 'Guiding wavelength = ', waveg
    return waveg
        

#====================================================

if __name__=="__main__":

    import sys
    import astropy.io.fits as fits
    import matplotlib.pyplot as plt
    import pylab as P
    
    inp = sys.argv[1]

    plt.close('all')

    #Set the temperature in K
    #Should be inherited from the main function in practice
    temp = 288.15

    # set zenith distance (degrees) and convert to airmass (= sec(z), with z in radians).
    #should also be inherited from main function
    z = 30.
    airmass = 1./np.cos(np.radians(z))

    # set humidity, pressure
    # Assume a relative humidity of 10\% (I used 9.6 for Paranal in the past)
    # Relative humidity = pw / ps * 100
    #(pw = partial pressure of water vapour, ps = saturation pressure of water vapour at the same temperature)
    #The value for ppwv below is valid for a RH = 10\%.
    #This and the atmosphere pressure of 760 mbar are values provided by ESO for Paranal.
    ppwv = 17.0
    p = 760.


    # read in the test cube
    hdu = fits.open(inp)
    cube = hdu[0].data
    head = hdu[0].header

    # from the FITS header keywords, create the wavelength axis. values are in MICRONS.
    wave = head['CRVAL3']+(np.arange(head['NAXIS3'])*head['CDELT3'])


    # create RA & dec axes also from the header, in MILLI-ARCSEC
    ra = np.arange(head['NAXIS1'])*head['CDELT1']
    dec = np.arange(head['NAXIS2'])*head['CDELT2']

    # pass to a function that calculates the differential refraction vector
    refr = calc_ref(wave, temp, ppwv, p)

    # now calculate the resulting differential refraction for the given airmass
    adr = calc_adr(wave, refr, airmass)

    #Apply ADR effect to datacube
    new_cube = add_ADR(cube, head, adr, 'linear')

    fits.writeto(sys.argv[1], new_cube, header=head)


##    P.figure(figsize=(11,9),dpi=72)
##    P.plot(wave, adr, 'k-', lw=1.5)
##    P.xlabel('$\lambda$ [$\mu$m]', fontsize=20)
##    P.ylabel('ADR [mas]', fontsize=20)
##    P.xticks(size=18)
##    P.yticks(size=18)
##    P.tick_params(which='major', length=12, width=1.5,
##                  axis='both')
##    P.tick_params(which='minor', length=8, width=1.5,
##                  axis='both')
##    P.show()
    
