'''Code that contains the functions dealing with
transmission and background aspects of the model

Written by Simon Zieleniewski

Created: 16-10-15

Last edited 01-12-15
'''

import numpy as n
import scipy as s
import os
from modules.gauss_convolve import gauss_convolve
from JSIM_Lowres import low_res_spec
from modules.misc_utils import path_setup


tppath = path_setup('../../Sim_data/Throughput/')



def create_thruput_cube(datacube_shape, wavels, resolution, grating, telescope=True, instrument=True,
                        illumination=None):
    '''Function that creates cube representing total throughput for each pixel in datacube.
    This incorporates sky, telescope, instrument and detector throughputs to convert photons to
    electrons. This will also be able to add in additional effects like illumination patterns.

    Inputs:
        datacube_shape: (z, y, x) shape of science datacube
        wavels: wavelength array of datacube
        resolution: spectral resolution [um]
        grating: Grating choice (string)
        telescope: Boolian incorporate telescope throughput
        instrument: Boolian - incorporate instrument throughput
        illumination: choice of illumination pattern


    Outputs:
        thruput_cube: cube representing total throughput of simulator. Science datacube
                      can be divided by this to give total throughput for each pixel.

    '''

    #Create emtpy cube of same size as science datacube
    init_cube = n.ones((datacube_shape), dtype=n.float64)

    if telescope:
        init_cube *= JWST_transmission_curve(wavels, grating)
        print 'Telescope transmission done!'
    if instrument:
        init_cube *= NIRSpec_transmission_curve(wavels)
        print 'Instrument transmission done!'


    #Add additional illumination and spectrograph patterns

    thruput_cube = init_cube

    print 'Throughput cube done!'

    return thruput_cube



#Telescope throughput curve generated just using wavelength array.
def JWST_transmission_curve(wavels, grating):
    '''Function that generates a full telescope throughput curve.

    Inputs:
        wavels: array of wavelengths for datacube
        grating: Grating choice (string)

    Outputs:
        cube_tele_trans: array of telescope throughput
                     for each wavelength value in wavels
    '''

    #Load grating throughput file
    #tele_r = n.genfromtxt(tppath+'EELT_mirror_reflectivity_mgf2agal.dat.txt')
    tele_r = n.genfromtxt(os.path.join(tppath,grating+'.csv'), delimiter=',')

    #Find start and end wavelength values of curves
    tt_start_arg = n.where(tele_r[:,0] < wavels[0])[0][-1]
    tt_end_arg = n.where(tele_r[:,0] > wavels[-1])[0][0]
    tele_trans_slice = tele_r[tt_start_arg:tt_end_arg+1,:]

    #Interpolate as a function of wavelength
    tele_trans_interp = s.interpolate.interp1d(tele_trans_slice[:,0], tele_trans_slice[:,1],
                                               kind='linear', bounds_error=False, fill_value=0.)
    #Obtain values for datacube wavelength array
    cube_tele_trans = tele_trans_interp(wavels)
    cube_tele_trans.shape = (len(wavels),1,1)

    return cube_tele_trans


#Instrument throughput curve generated just using wavelength array.
def NIRSpec_transmission_curve(wavels):
    '''Function that generates a full NIRSpec IFU throughput curve

    Inputs:
        wavels: array of wavelengths for datacube

    Outputs:
        cube_inst_trans: array of instrument throughput
                     for each wavelength value in wavels
    '''

    #Load telescope reflectivity file
    inst_r = n.genfromtxt(os.path.join(tppath,'NIRSpecthruput.txt'), delimiter=',')

    #Find start and end wavelength values of curves
    it_start_arg = n.where(inst_r[:,0] < wavels[0])[0][-1]
    it_end_arg = n.where(inst_r[:,0] > wavels[-1])[0][0]
    inst_trans_slice = inst_r[it_start_arg:it_end_arg+1,:]

    #Interpolate as a function of wavelength
    inst_trans_interp = s.interpolate.interp1d(inst_trans_slice[:,0], inst_trans_slice[:,1],
                                               kind='linear', bounds_error=False, fill_value=0.)
    #Obtain values for datacube wavelength array
    cube_inst_trans = inst_trans_interp(wavels)
    cube_inst_trans.shape = (len(wavels),1,1)

    return cube_inst_trans
