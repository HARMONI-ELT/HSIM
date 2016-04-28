'''Code that contains the functions dealing with
transmission and background aspects of the model

Written by Simon Zieleniewski

Started 11-06-13
Last edited 27-04-16
'''

import numpy as n
import scipy as s
import os
from modules.gauss_convolve import gauss_convolve
from TheSimulator_Lowres import low_res_spec
from modules.misc_utils import path_setup


tppath = path_setup('../../Sim_data/Throughput/')



def create_thruput_cube(datacube_shape, wavels, resolution, grating, inst_tpvals, sky=True, telescope=True, instrument=True, QE=True,
                        illumination=None):
    '''Function that creates cube representing total throughput for each pixel in datacube.
    This incorporates sky, telescope, instrument and detector throughputs to convert photons to
    electrons. This will also be able to add in additional effects like illumination patterns.

    Inputs:
        datacube_shape: (z, y, x) shape of science datacube
        wavels: wavelength array of datacube
        resolution: spectral resolution [um]
        grating: grating choice to set grating throughput curve
        inst_tpvals: instrument throughput values (from config file) [W_grat, WO_grat]
        sky: Boolian - incorporate sky throughput
        telescope: Boolian incorporate telescope throughput
        instrument: Boolian - incorporate instrument throughput
        QE: Boolian - incorporate detector QE throughput
        illumination: choice of illumination pattern


    Outputs:
        thruput_cube: cube representing total throughput of simulator. Science datacube
                      can be divided by this to give total throughput for each pixel.

    '''

    #Create emtpy cube of same size as science datacube
    init_cube = n.ones((datacube_shape), dtype=n.float64)

    if sky:
        init_cube *= sky_transmission_curve(wavels, resolution)
        print 'Sky tranismission done!'
    if telescope:
        init_cube *= telescope_transmission_curve(wavels)
        print 'Telescope transmission done!'
    if instrument:
        init_cube *= HARMONI_transmission_curve(wavels, grating, inst_tpvals)
        print 'Instrument transmission done!'
    if QE:
        init_cube *= detector_QE_curve(wavels)
        print 'Detector QE done!'

    #Add additional illumination and spectrograph patterns

    thruput_cube = init_cube

    print 'Throughput cube done!'

    return thruput_cube




#Sky throughput curve generated just using wavelength array.
def sky_transmission_curve(wavels, delta_lambda):
    '''Function that generates a full throughput curve combining
    sky transmission & sky extinction.

    Inputs:
        wavels: array of wavelengths for datacube
        delta_lambda: Resolution element [um]

    Outputs:
        cube_total_sky_trans: array of total throughput
                     for each wavelength value in wavels
    '''

    #NEED TO IMPLEMENT DIFFERENT AIRMASS VALUES AT SOME POINT
    #Currently set at 1.15 to match PSFs - Zenith angle = 30deg

    #load sky transmission & extinction files
    sky_trans = n.genfromtxt(os.path.join(tppath,'ASM_throughput/transmission_resolution_0.15_angstroms_MICRONS.txt'))

    #Find start and end wavelength values of curves
    st_start_arg = n.where(sky_trans[:,0] < wavels[0])[0][-1]
    st_end_arg = n.where(sky_trans[:,0] > wavels[-1])[0][0]
    sky_trans_slice = sky_trans[st_start_arg:st_end_arg+1,:]

    if type(delta_lambda) == n.ndarray:
        #Low resolution mode. Use low_res_mode function to generate sky transmission spectrum
        sky_trans_spec = low_res_spec(sky_trans, delta_lambda, True)
        sky_total_trans = sky_trans_spec[:,1].copy()

    else:
        #Convolve sky transmission array with Gaussian LSF of
        #FWHM = sqrt(new_resolution**2-old_resolution**2)
        #to fold in spectrum for each spectral pixel.
        sigma = n.sqrt(delta_lambda**2-0.15E-4**2)/(2.*n.sqrt(2.*n.log(2.)))
        conv_sky_trans = gauss_convolve(sky_trans_slice, sigma, lambda_space='Linear')

        interp_trans = s.interpolate.interp1d(conv_sky_trans[:,0], conv_sky_trans[:,1],
                                              kind='linear')
        sky_total_trans = interp_trans(wavels)

    sky_total_trans.shape = (len(wavels),1,1)
    return sky_total_trans


#Telescope throughput curve generated just using wavelength array.
def telescope_transmission_curve(wavels):
    '''Function that generates a full telescope throughput curve.
    Current telescope design contains 6 mirros, each with the same
    Ag/Al coating and reflectivity curve. Thus total throughput
    given by:
    T(lambda) = R(lambda)^6

    Inputs:
        wavels: array of wavelengths for datacube

    Outputs:
        cube_tele_trans: array of telescope throughput
                     for each wavelength value in wavels
    '''

    #Telescope transmission model taken from EELT DRM telescope model
    #http://www.eso.org/sci/facilities/eelt/science/drm/tech_data/telescope/

    #Load telescope reflectivity file
    tele_r = n.genfromtxt(os.path.join(tppath,'EELT_mirror_reflectivity_mgf2agal.dat.txt'))

    #Apply dust free area fraction as calculated by Niranjan Thatte 02-06-14
    tele_r[:,1] *= 0.978

    #Find start and end wavelength values of curves
    tt_start_arg = n.where(tele_r[:,0] < wavels[0])[0][-1]
    tt_end_arg = n.where(tele_r[:,0] > wavels[-1])[0][0]
    tele_trans_slice = tele_r[tt_start_arg:tt_end_arg+1,:]

    #Interpolate as a function of wavelength
    tele_trans_interp = s.interpolate.interp1d(tele_trans_slice[:,0], tele_trans_slice[:,1],
                                               kind='linear', bounds_error=False, fill_value=0.)
    #Obtain values for datacube wavelength array
    cube_tele_ref = tele_trans_interp(wavels)

    #Apply equation to obtain transmission (7 mirror design)
    cube_tele_trans = cube_tele_ref ** 7.
    cube_tele_trans.shape = (len(wavels),1,1)

    return cube_tele_trans


#Instrument throughput curve generated just using wavelength array.
def HARMONI_transmission_curve(wavels, grating, inst_tpvals):
    '''Function that generates a full HARMONI throughput curve.
    Combines grating throughput curve with flat instrument value

    Inputs:
        wavels: array of wavelengths for datacube
        grating: grating choice to set grating throughput curve
        inst_tpvals: instrument throughput values (from config file) [W_grat, WO_grat]

    Outputs:
        cube_inst_trans: array of instrument throughput
                     for each wavelength value in wavels
    '''

    # #Load instrument throughput file
    # inst_r = n.genfromtxt(os.path.join(tppath,'HARMONIthruput.txt'), delimiter=',')
    #
    # #Find start and end wavelength values of curves
    # it_start_arg = n.where(inst_r[:,0] < wavels[0])[0][-1]
    # it_end_arg = n.where(inst_r[:,0] > wavels[-1])[0][0]
    # inst_trans_slice = inst_r[it_start_arg:it_end_arg+1,:]
    #
    # #Interpolate as a function of wavelength
    # inst_trans_interp = s.interpolate.interp1d(inst_trans_slice[:,0], inst_trans_slice[:,1],
    #                                            kind='linear', bounds_error=False, fill_value=0.)
    # #Obtain values for datacube wavelength array
    # cube_inst_trans = inst_trans_interp(wavels)
    # cube_inst_trans.shape = (len(wavels),1,1)

    if grating != 'None':
        if grating in ['V+R', 'V', 'R', 'V-high', 'R-high']:
            gratingfile = 'V+R'
        elif grating in ['Iz+J', 'Iz', 'J', 'z', 'J-high']:
            gratingfile = 'Iz+J'
        elif grating in ['H+K', 'H', 'K', 'H-high', 'K-high']:
            gratingfile = 'H+K'
        #Load grating throughput file
        inst_r = n.genfromtxt(os.path.join(tppath,gratingfile+'_grating.txt'), delimiter=',')

        #Find start and end wavelength values of curves
        it_start_arg = n.where(inst_r[:,0] < wavels[0])[0][-1]
        it_end_arg = n.where(inst_r[:,0] > wavels[-1])[0][0]
        inst_trans_slice = inst_r[it_start_arg:it_end_arg+1,:]

        #Interpolate as a function of wavelength
        inst_trans_interp = s.interpolate.interp1d(inst_trans_slice[:,0], inst_trans_slice[:,1],
                                                   kind='linear', bounds_error=False, fill_value=0.)
        #Obtain values for datacube wavelength array
        cube_inst_trans = inst_trans_interp(wavels)*inst_tpvals[0]
        cube_inst_trans.shape = (len(wavels),1,1)

    elif grating == 'None':
        cube_inst_trans = n.ones(len(wavels), dtype=float)*inst_tpvals[1]
        cube_inst_trans.shape = (len(wavels),1,1)

    return cube_inst_trans


#Detector throughput curve generated just using wavelength array.
def detector_QE_curve(wavels):
    '''Function that generates a detector QE curve.
    Current detector data taken from Detector Reference
    table.

    Inputs:
        wavels: array of wavelengths for datacube

    Outputs:
        cube_det_qe: array of detector QE values for
                  each wavelength in array
    '''

    #Load telescope reflectivity file
    det_qe = n.genfromtxt(os.path.join(tppath,'DetectorQE.txt'), delimiter=',')

    #Find start and end wavelength values of curves
    qe_start_arg = n.where(det_qe[:,0] < wavels[0])[0][-1]
    qe_end_arg = n.where(det_qe[:,0] > wavels[-1])[0][0]
    det_qe_slice = det_qe[qe_start_arg:qe_end_arg+1,:]

    #Interpolate as a function of wavelength
    det_qe_interp = s.interpolate.interp1d(det_qe_slice[:,0], det_qe_slice[:,1],
                                           kind='linear', bounds_error=False, fill_value=0.)
    #Obtain values for datacube wavelength array
    cube_det_qe = det_qe_interp(wavels)
    cube_det_qe.shape = (len(wavels),1,1)

    return cube_det_qe
