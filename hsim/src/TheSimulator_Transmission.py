'''Code that contains the functions dealing with
transmission and background aspects of the model

Written by Simon Zieleniewski

Started 11-06-13
Last edited 27-10-16 by Laurence Routledge
'''

import numpy as n
import scipy as s
import os
from modules.gauss_convolve import gauss_convolve
from TheSimulator_Lowres import low_res_spec
from modules.misc_utils import path_setup


tppath = path_setup('../../Sim_data/Throughput/')

def create_thruput_cube(datacube_shape, wavels, resolution, grating, zenith_angle, inst_tpvals, sky=True, telescope=True,
			instrument=True, QE=True, illumination=None):
	'''Function that creates cube representing total throughput for each pixel in datacube.
	This incorporates sky, telescope, instrument and detector throughputs to convert photons to
	electrons. This will also be able to add in additional effects like illumination patterns.

	Inputs:
		datacube_shape: (z, y, x) shape of science datacube
		wavels: wavelength array of datacube
		resolution: spectral resolution [um]
		grating: grating choice to set grating throughput curve
		zenith_angle: zenith angle of observation [degrees]
		inst_tpvals: instrument throughput values (from config file) [WO_grat]
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
		sky_cube, wavels = sky_transmission_curve(wavels, resolution, zenith_angle)
		init_cube *= sky_cube
		print 'Sky tranismission done!'
	if telescope:
		init_cube *= telescope_transmission_curve(wavels)
		print 'Telescope transmission done!'
	if instrument:
		init_cube *= HARMONI_transmission_curve(wavels, grating, inst_tpvals)
		print 'Instrument transmission done!'
	if QE:
		init_cube *= detector_QE_curve(wavels, grating)
		print 'Detector QE done!'

	#Add additional illumination and spectrograph patterns

	thruput_cube = init_cube

	print 'Throughput cube done!'

	return thruput_cube




#Sky throughput curve generated just using wavelength array.
def sky_transmission_curve(wavels, delta_lambda, zenith_ang):
	'''Function that generates a full throughput curve combining
	sky transmission & sky extinction.

	Inputs:
		wavels: array of wavelengths for datacube
		delta_lambda: Resolution element [um]
		zenith_ang: zenith angle of observation [degrees]

	Outputs:
		cube_total_sky_trans: array of total throughput
			for each wavelength value in wavels
	'''
	#convert from zenith angle to airmass
	X = 1./n.cos(n.radians(zenith_ang))
	inbuilt_airmasses = [1., 1.15, 1.41, 2.]

	#determine the closest data to the airmass value given and find it's location in the data file
	closest_X = min(inbuilt_airmasses, key=lambda x:abs(x - X))
	data_index = inbuilt_airmasses.index(closest_X) + 1

	#load sky transmission & extinction files, then reduce to the columns required
	sky_trans_all_X = n.genfromtxt(os.path.join(tppath,'ASM_throughput/transmission_0.15_angstroms_resolution.txt'))
	sky_trans = n.column_stack((sky_trans_all_X[:,0], sky_trans_all_X[:,data_index]))
	print len(sky_trans)
	sky_trans_all_X = []            #clean up unused data for sake of memory
		
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
				kind='linear', bounds_error=False, fill_value=0)
		sky_total_trans = interp_trans(wavels)

	sky_total_trans.shape = (len(wavels),1,1)
	return sky_total_trans, wavels


#Telescope throughput curve generated just using wavelength array.
def telescope_transmission_curve(wavels):
	'''Function that reads a full telescope throughput curve.

Inputs:
	wavels: array of wavelengths for datacube

Outputs:
	cube_tele_trans: array of telescope throughput
		for each wavelength value in wavels
	'''

	#Load telescope reflectivity file
	tele_r = n.genfromtxt(os.path.join(tppath,'ELT_mirror_reflectivity.txt'), delimiter=',')

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
def HARMONI_transmission_curve(wavels, grating, inst_tpvals):
	'''Function that generates a full HARMONI throughput curve.
	Combines grating throughput curve with flat instrument value

	Inputs:
		wavels: array of wavelengths for datacube
		grating: grating choice to set grating throughput curve
		inst_tpvals: instrument throughput values (from config file) [ WO_grat]

	Outputs:
		cube_inst_trans: array of instrument throughput
			for each wavelength value in wavels
	'''

	if grating != 'None':
		gratingfile = grating

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
		cube_inst_trans = inst_trans_interp(wavels)
		cube_inst_trans.shape = (len(wavels),1,1)

	elif grating == 'None':
		cube_inst_trans = n.ones(len(wavels), dtype=float)*inst_tpvals[0]
		cube_inst_trans.shape = (len(wavels),1,1)

	return cube_inst_trans


#Detector throughput curve generated just using wavelength array.
def detector_QE_curve(wavels, grating):
	'''Function that generates a detector QE curve.
	Current detector data taken from Detector Reference
	table.

	Inputs:
		wavels: array of wavelengths for datacube
		grating: grating choice to set detector

	Outputs:
		cube_det_qe: array of detector QE values for
			each wavelength in array
	'''

	if grating == "V+R":
		detector_QE_file = "Psyche_CCD231-84_ESO_measured_QE.txt"
	else:
		detector_QE_file = "H4RG_QE_design.txt"
		
	det_qe = n.genfromtxt(os.path.join(tppath,detector_QE_file), delimiter=',')

	#Find start and end wavelength values of curves
	qe_start_arg = n.where(det_qe[:,0] < wavels[0])[0][-1]
	qe_end_arg = n.where(det_qe[:,0] > wavels[-1])[0][0]
	det_qe_slice = det_qe[qe_start_arg:qe_end_arg+1,:]

	#Interpolate as a function of wavelength
	det_qe_interp = s.interpolate.interp1d(det_qe_slice[:,0], det_qe_slice[:,1],
				kind='linear', bounds_error=False, fill_value=0.)
	#Obtain values for datacube wavelength array
	cube_det_qe = det_qe_interp(wavels)/100.
	cube_det_qe.shape = (len(wavels),1,1)

	return cube_det_qe
