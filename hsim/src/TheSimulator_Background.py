'''Code that contains the functions dealing with
the sky & thermal background aspects of the model.

Written by Simon Zieleniewski

Started 14-06-13
Last edited: 23-06-15
'''


import os
import numpy as n
from scipy.interpolate import interp1d
import scipy.constants as sc
from modules.gauss_convolve import gauss_convolve
from modules.frebin import *
from TheSimulator_Lowres import low_res_spec
from modules.misc_utils import path_setup


bgpath = path_setup('../../Sim_data/Background/')


def create_background_cube(datacube_shape, wavels, thru_cube, inst_qe_cube, qe_cube, DIT, NDIT, spaxel,
                           resolution, wave_disp, area, temp, sky=True, telescope=True, instrument=True,
                           emis=0.244, slow=False):
    '''Function that creates cube representing total throughput for each pixel in datacube.
    This incorporates sky, telescope, instrument and detector throughputs to convert photons to
    electrons. This will also be able to add in additional effects like illumination patterns and
    separate spectrograph quartiles throughput.

    Inputs:
        datacube_shape: (z, y, x) shape of science datacube
        wavels: wavelength array of datacube
        thru_cube: cube representing total throughput.
                   Used to convert sky signal from photons > electrons
        tel_qe_cube: cube representing throughput of instrument+detector.
                     Used to convert telescope thermal background from photons > electrons
        qe_cube: cube representing quantum effiency of detector.
                 Used to convert signal from telescope and instrument from photons > electrons.
        DIT: exposure time [s]
        NDIT: no. of exposures
        spaxel: spaxel scale (x, y) [mas, mas]
        resolution: resolution element [um]
        wave_disp: wavelength sampling [um/pix]
        area: telescope area [m]
        temp: Telescope temperature [K]
        sky: Boolian - incorporate sky continuum, emission lines and OH lines
        telescope: Boolian incorporate telescope thermal background
        instrument: Boolian - incorporate instrument thermal background
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
    if instrument:
        init_cube += HARMONI_background(wavels)*qe_cube
        print 'Instrument background done!'
        

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
    '''Function that generates a sky background curve combining
    sky continuum, sky thermal emission and sky emission lines.
    
    Inputs:
        wavels: array of wavelengths for datacube
        delta_lambda: resolution element [um]
        dit: exposure time [s]. This determins how the sky emission
             line amplitudes vary through the exposure.
        
    Outputs:
        sky_bg_curve: array of total sky background
                  [units of photons/s/m^2/um/arcsec^2]
                  for each wavelength value in wavels
    '''
    
    #Load up all datafile using n.genfromtxt(...)
    #sky_em = n.genfromtxt(bgpath+'ASM_background/radiance_resolution_0.15_angstroms_MICRONS.txt')
    sky_em = n.genfromtxt(os.path.join(bgpath,'ASM_background/radiance_resolution_0.15_angstroms_MICRONS.txt'))
    
    #Cut down data to relevent region.
    se_start_arg = n.where(sky_em[:,0] < wavels[0])[0][-1]
    se_end_arg = n.where(sky_em[:,0] > wavels[-1])[0][0]
    sky_em_slice = sky_em[se_start_arg:se_end_arg+1,:]

    if type(delta_lambda) == n.ndarray:
        #Low resolution mode. Use low_res_mode function to generate sky emission spectrum
        sky_em_spec = low_res_spec(sky_em, delta_lambda, False)
        binned_sky_em = sky_em_spec[:,1].copy()

    else:
        #Convolve sky transmission array with Gaussian LSF of
        #FWHM = sqrt(new_resolution**2-old_resolution**2)
        #to fold in spectrum for each spectral pixel.
        input_disp = sky_em[1,0] - sky_em[0,0]
        sigma = n.sqrt(delta_lambda**2-input_disp**2)/(2.*n.sqrt(2.*n.log(2.)))
        
        conv_sky_em = gauss_convolve(sky_em_slice, sigma, lambda_space='Linear')
        ###new 30-05-14
        conv_sky_em_ph = conv_sky_em[:,1]*input_disp #convert to units of photons/s/m^2/arcsec^2
        binned_sky_em_ph = frebin(conv_sky_em_ph.reshape(len(conv_sky_em_ph),1), (1,len(wavels)), True) #Total photons conserved
        binned_sky_em = binned_sky_em_ph/float(wavels[1]-wavels[0]) #reconvert back to photons/s/m^2/um/arcsec^2
    
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
    #The EELT will have (at least) 6 mirros before entering the
    #HARMONI cryostat (at which point thermal background is negligible).
    #Average Paranal night temperature, T = 285K (from ESO website).
    
    #Blackbody function
    cube_bb_spec = blackbody(wavels, T)

    tele_bg_spec = emissivity*cube_bb_spec
    tele_bg_spec_ph = tele_bg_spec/(sc.h*sc.c/(wavels*1.E-6))
    tele_bg_spec_ph.shape = (len(wavels),1,1)
    return tele_bg_spec_ph


def HARMONI_background(wavels):
    '''Function that generates a HARMONI background curve
    using mirror reflectivities, mirror areas and the
    Plank BB function as a function of instrument temperature.
    HARMONI modelled as a blackbody BB(T) with T = 130 K.
    
    Inputs:
        wavels: array of wavelengths for datacube
        D: telescope diameter [m]
        area: telescope area [m^2]
        eps: telescope obscuration ratio
        
    Outputs:
        tele_bg_curve: array of total telescope background
                  [units of photons/s/m^2/um/arcsec^2]
                  for each wavelength value in wavels
    '''

    #Instrument at 130 K
    #Blackbody function
    cube_bb_spec = blackbody(wavels, 130.0)

    inst_bg_spec = cube_bb_spec
    inst_bg_spec_ph = inst_bg_spec/(sc.h*sc.c/(wavels*1.E-6))
    inst_bg_spec_ph.shape = (len(wavels),1,1)
    return inst_bg_spec_ph
