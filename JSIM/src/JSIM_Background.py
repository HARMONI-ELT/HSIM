'''Code that contains the functions dealing with
the sky & thermal background aspects of the model.

Written by Simon Zieleniewski

Started 16-10-15
Last edited: 01-12-15
'''


import os
import numpy as n
from scipy.interpolate import interp1d
import scipy.constants as sc
from modules.gauss_convolve import gauss_convolve
from modules.frebin import *
from JSIM_Lowres import low_res_spec
from modules.misc_utils import path_setup


bgpath = path_setup('../../Sim_data/Background/')


def create_background_cube(datacube_shape, wavels, thru_cube, inst_qe_cube, DIT, NDIT, spaxel,
                           resolution, wave_disp, area, temp, sky=True, telescope=True,
                           emis=0.15, slow=False):
    '''Function that creates cube representing total background photons for each pixel in datacube.

    Inputs:
        datacube_shape: (z, y, x) shape of science datacube
        wavels: wavelength array of datacube
        thru_cube: cube representing total throughput.
                   Used to convert sky signal from photons > electrons
        tel_qe_cube: cube representing throughput of instrument+detector.
                     Used to convert telescope thermal background from photons > electrons
        DIT: exposure time [s]
        NDIT: no. of exposures
        spaxel: spaxel scale (x, y) [mas, mas]
        resolution: resolution element [um]
        wave_disp: wavelength sampling [um/pix]
        area: telescope area [m]
        temp: Telescope temperature [K]
        sky: Boolian - incorporate sky continuum, emission lines and OH lines
        telescope: Boolian incorporate telescope thermal background
        emis: telescope emissivity - default is 24.4% after 7 warm reflections and with unmasked spider arms.


    Outputs:
        background_cube_new: cube representing total background of simulator.
        background_noise: cube representing shot noise on background
        background_cube: cube representing noiseless total background

    '''

    print 'Generating background cube.'

    #Create emtpy cube of same size as science datacube
    init_cube = n.zeros((datacube_shape), dtype=n.float64)

    if sky:
        init_cube += sky_background(wavels, resolution, DIT)*thru_cube
        print 'Sky background done!'
    if telescope:
        init_cube += telescope_background(wavels, temp, emis)*inst_qe_cube
        print 'Telescope background done!'


    #Add additional effects here

    background_cube = init_cube

    if slow==True:
        background_cube *= DIT*(spaxel[0]*spaxel[1])*wave_disp*area

        background_cube_new = n.zeros(background_cube.shape, dtype=float)
        for x in xrange(NDIT):
            #Shot noise
            background_cube_new += n.random.poisson(background_cube)
        background_cube_new/=n.float(NDIT)

    else:
        #Generate total background cube
        background_cube *= DIT*NDIT*(spaxel[0]*spaxel[1])*wave_disp*area
        #Generate background noise
        background_cube_new = n.random.poisson(n.abs(background_cube))

    background_cube_new = background_cube_new.astype(n.float64)
    background_noise = n.abs(background_cube_new-background_cube)

    print 'Background cube generated!'

    return background_cube_new, background_noise, background_cube



def blackbody(waves, T):
    '''Function that gives the Plank blackbody spectrum
    as a function of wavelength and temperature.

    Inputs:
        waves: value or array of wavelengths in microns
        T: Temperature in Kelvin

    Outputs:
        bb_spectrum: value or array for Plank blackbody spectrum
                     units - [J/s/m^2/lambda(um)/arcsec^2]
    '''

    #Convert wavelength array from microns to metres
    wave = waves*1.E-6

    exp_part = n.exp(sc.h*sc.c/(wave*sc.k*T))

    #Flux in [J/s/m/m2/steradian]
    bb_spectrum = (2.*sc.h*sc.c**2/wave**5)*(exp_part - 1)**(-1)
    #put into units of: J/s/lambda(um)/m^2/arcsec^2
    bb_spectrum /= 1.E6 #to get into J/s/m2/um/steradian
    bb_spectrum /= 4.2545E10 #to get into J/s/m2/um/arcsec2


    return bb_spectrum




def sky_background(wavels, delta_lambda, dit):
    '''Function that generates a sky background curve using
    Zodiacal light. From JWST ETC.

    Inputs:
        wavels: array of wavelengths for datacube
        dit: exposure time [s]. This determins how the sky emission
             line amplitudes vary through the exposure.

    Outputs:
        sky_bg_curve: array of total sky background
                  [units of photons/s/m^2/um/arcsec^2]
                  for each wavelength value in wavels
    '''

    #Load up all datafile using n.genfromtxt(...)
    #File in units [J/s/m2/um/arcsec2]
    sky_em = n.genfromtxt(os.path.join(bgpath,'Zodiacal_background_cgs.csv'), delimiter=',')

    #Cut down data to relevent region.
    se_start_arg = n.where(sky_em[:,0] < wavels[0])[0][-1]
    se_end_arg = n.where(sky_em[:,0] > wavels[-1])[0][0]
    sky_em_slice = sky_em[se_start_arg:se_end_arg+1,:]

    if type(delta_lambda) == n.ndarray:
        #Low resolution mode. Use low_res_mode function to generate sky emission spectrum
        sky_em_spec = low_res_spec(sky_em, delta_lambda, False)
        binned_sky_em = sky_em_spec[:,1].copy()

    else:
        #Interpolate as a function of wavelength
        sky_em_interp = interp1d(sky_em_slice[:,0], sky_em_slice[:,1],
                                                   kind='linear', bounds_error=False, fill_value=0.)
        #Obtain values for datacube wavelength array
        binned_sky_em = sky_em_interp(wavels)
        binned_sky_em.shape = (len(wavels),1,1)

    return binned_sky_em


def telescope_background(wavels, T, emissivity):
    '''Function that generates a telescope background curve
    using mirror reflectivities, mirror areas and the
    Plank BB function as a function of site temperature.
    Telescope modelled as a graybody with eps*BB(T), where
    eps = emissivity of telescope. Emissivity is computed with
    mirror reflectivites for 6 mirrors of telescope.

    Inputs:
        wavels: array of wavelengths for datacube
        T: site temperature [K]

    Outputs:
        tele_bg_curve: array of total telescope background
                  [units of photons/s/m^2/um/arcsec^2]
                  for each wavelength value in wavels
    '''

    #TELESCOPE emission should be modelled as a graybody: blackbody
    #multiplied by a constant emissivity. Emissivity is calculated
    #by (1-R)^(no. of mirros) where R is the mirror reflectivity.

    #Blackbody function
    cube_bb_spec = blackbody(wavels, T)

    tele_bg_spec = emissivity*cube_bb_spec
    tele_bg_spec_ph = tele_bg_spec/(sc.h*sc.c/(wavels*1.E-6))
    tele_bg_spec_ph.shape = (len(wavels),1,1)
    return tele_bg_spec_ph
