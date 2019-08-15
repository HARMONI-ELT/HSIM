'''
Calculates detector QE, dark current, read noise
'''
import os
import logging

import numpy as np

import scipy.constants as sp
from scipy.interpolate import UnivariateSpline
from scipy import integrate
from astropy.convolution import Gaussian1DKernel

from src.modules.misc_utils import path_setup
from src.config import *
from src.modules.em_model import *

detpath = path_setup('../../' + config_data["data_dir"] + 'detectors/')

#Detector throughput curve generated just using wavelength array.
def detector_QE_curve(wavels, grating, debug_plots, output_file):
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
		
		
	cube_det_qe, orig_qe_lambda, orig_qe = load_transmission_curve(wavels, detector_QE_file, debug_plots, [output_file, "det_qe"], "detector QE", scaling=1/100., full_curve=True)
	return cube_det_qe, orig_qe_lambda, orig_qe

def mask_saturated_pixels(cube, grating):
	''' Mask saturated pixels
	Inputs:
		cube: Input datacube (RA, DEC, lambda)
		grating: Spectral grating
	Outputs:
		cube: masked cube
		mask_pixels: saturated pixels mask
	'''
	logging.info("Masking saturated pixels")
	
	if grating == "V+R":
		limit = config_data["saturation"]["vis"]
	else:
		limit = config_data["saturation"]["nir"]
		
	
	mask_pixels = cube > limit
	if np.sum(mask_pixels) > 0:
		logging.warning(str(np.sum(mask_pixels)) + " pixels are saturated in the observed cube")
		
	cube[mask_pixels] = 0.
	
	return cube, mask_pixels

	

def apply_crosstalk(cube, crosstalk):
	''' Simulates crosstalk detector effects
	Inputs:
		cube: Input datacube (RA, DEC, lambda)
		crosstalk: Fraction of the photons that go to each of the 4 contiguous pixels
	Outputs:
		cube: Cube including crosstalk
	'''
	
	#logging.info("Applying detector crosstalk")
	
	scaled_cube = cube*(1. - crosstalk*4)
	
	# crosstalk in the spatial direction
	spatial_crosstalk = crosstalk*(np.roll(cube, 1, axis=2) + np.roll(cube, -1, axis=2))
	
	# crosstalk in the spectral direction
	spectral_crosstalk = crosstalk*(np.roll(cube, 1, axis=0) + np.roll(cube, -1, axis=0))
	
	return scaled_cube + spatial_crosstalk + spectral_crosstalk


def apply_crosstalk_1d(spectrum, crosstalk):
	''' Simulates crosstalk detector effects
	Inputs:
		spectrum: Input spectrum (lambda)
		crosstalk: Fraction of the photons that go to each of the 4 contiguous pixels
	Outputs:
		spectrum including crosstalk
	'''
	
	#logging.info("Applying detector crosstalk - 1d")
	
	scaled_spectrum = spectrum*(1. - crosstalk*2)
	
	# crosstalk in the spectral direction
	spectral_crosstalk = crosstalk*(np.roll(spectrum, 1, axis=0) + np.roll(spectrum, -1, axis=0))
	
	return scaled_spectrum + spectral_crosstalk


 

def sim_detector(cube, back_emission, transmission, lambs, grating, DIT, debug_plots=False, output_file=""):
	''' Simulates detector effects
	Inputs:
		cube: Input datacube (RA, DEC, lambda)
		back_emission: Input background emission
		transmission: Input transmission
		lambs: lambda array [um]
		grating: Spectral grating
		debug_plots: Produce debug plots
		output_file: File name for debug plots

	Outputs:
		cube: Cube including QE
		read_noise: read noise for the grating and DIT [e/pix]
		dark*DIT: dark current  [e/pix]
		ftot_electron*DIT: thermal background seen by the detector  [e/pix]
	'''
	
	# Get QE curve
	logging.info("Calculating detector QE")
	qe_curve, orig_qe_lambda, orig_qe = detector_QE_curve(lambs, grating, debug_plots, output_file)
	back_emission *= qe_curve
	transmission *= qe_curve
	
	qe_curve.shape = (len(lambs), 1, 1)
	cube *= qe_curve

	# read noise
	if grating == "V+R":
		read_noise = config_data["read_noise"]["vis"]
	else:
		if DIT <= 120.:
			read_noise = config_data["read_noise"]["nir_lowexp"]
		else:
			read_noise = config_data["read_noise"]["nir"]


	# dark current
	if grating == "V+R":
		dark = config_data["dark_current"]["vis"]
	else:
		dark = config_data["dark_current"]["nir"]


	# Thermal emission seen by the detector
	# 15-micron pixels and F/2 camera... (independent of spaxel scale)
	pixel_area = 15E-6**2 # m2
	pixel_solid_cryo_mech = 0.2*(360/2./np.pi*3600.)**2 # arcsec**2
	pixel_solid_rad = 2.0*np.pi*(360/2./np.pi*3600.)**2 # arcsec**2
	
	TCryoMech = config_data["HARMONI_cryo_temp"] - 5.
	Trad = config_data["HARMONI_cryo_temp"]
	
	fcryomech = blackbody(orig_qe_lambda, TCryoMech)/(sp.h*sp.c/(orig_qe_lambda*1.E-6))*pixel_solid_cryo_mech # photons/s/um/m2
	frad = blackbody(orig_qe_lambda, Trad)/(sp.h*sp.c/(orig_qe_lambda*1.E-6))*pixel_solid_rad # photons/s/um/m2
	ftot = (fcryomech + frad)*pixel_area # photons/um/pixel
	ftot = ftot*orig_qe  # e/um/pixel
	ftot_electron = integrate.simps(ftot, orig_qe_lambda) #e/pixel
	
	return cube, back_emission, transmission, read_noise, dark*DIT, ftot_electron*DIT
	

# Detector systematics code

#interpolation function
def interp(data):
	X = np.sort(data)
	x_range = np.linspace(0,1,len(X))
	func = UnivariateSpline(x_range, X, k=4, s=0)
	return func


def make_det_vals(data, N):
	np.random.seed(1)
	rands = np.random.random(N)
	rands_sort = np.sort(rands)
	det_vals = interp(data)(rands_sort)
	np.random.seed(1)
	np.random.shuffle(det_vals)
	return det_vals


def make_detectors(NDIT):
	''' Generates an instance of the detectors
	Inputs:
		NDIT: Number of detector integrations

	Outputs:
		sim_dets: list of simulated detectors
	'''
	# load KMOS maps
	means = np.load(os.path.join(detpath, 'kmos_means.npy'))
	var = np.load(os.path.join(detpath, 'kmos_var.npy'))
	side_len = config_data["side_length"]

	sds = np.sqrt(var)

	mc = means[4:-4,4:-4].flatten()
	mt = means[:,-4:].flatten()
	ml = means[4:-4,:4].flatten() 
	mr = means[4:-4,-4:].flatten() 
	
	sc = sds[4:-4,4:-4].flatten()
	sc[sc > 1000] = 1000
	st = sds[:,-4:].flatten()
	sl = sds[4:-4,:4].flatten()
	sr = sds[4:-4,-4:].flatten()
	
	mean_map = np.zeros((8,side_len,side_len))
	sds_map = np.zeros((8,side_len,side_len))

	mean_map[:,4:-4,4:-4] = make_det_vals(mc, ((side_len-8)**2)*8).reshape((8,side_len-8, side_len-8))
	mean_map[:,-4:,:] = make_det_vals(mt, (side_len*4)*8).reshape((8,4,side_len))
	mean_map[:,:4,:] = make_det_vals(mt, (side_len*4)*8).reshape((8,4,side_len))
	mean_map[:,4:-4,:4] = make_det_vals(ml, ((side_len-8)*4)*8).reshape((8,side_len-8,4))
	mean_map[:,4:-4,-4:] = make_det_vals(mr, ((side_len-8)*4)*8).reshape((8,side_len-8,4))

	sds_map[:,4:-4,4:-4] = make_det_vals(sc, ((side_len-8)**2)*8).reshape((8,side_len-8, side_len-8))
	sds_map[:,-4:,:] = make_det_vals(st, (side_len*4)*8).reshape((8,4,side_len))
	sds_map[:,:4,:] = make_det_vals(st, (side_len*4)*8).reshape((8,4,side_len))
	sds_map[:,4:-4,:4] = make_det_vals(sl, ((side_len-8)*4)*8).reshape((8,side_len-8,4))
	sds_map[:,4:-4,-4:] = make_det_vals(sr, ((side_len-8)*4)*8).reshape((8,side_len-8,4))
	return mean_map, sds_map


def make_det_instance(NDIT):
	mean_map, sd_map = make_detectors(NDIT)
	side_len = config_data["side_length"]
	sim_dets = sd_map * np.random.randn(8,side_len,side_len) + mean_map
	col_meds = np.load(os.path.join(detpath, 'col_meds.npy'))

	for i in range(mean_map.shape[0]):
		for j in range(mean_map.shape[1]-8):
			row=j+4
			row_med = np.median(np.concatenate([sim_dets[i,row-1:row+2,:4], sim_dets[i,row-1:row+2,-4:]]))
			sim_dets[i,row] -= row_med

		for k in range(mean_map.shape[2]/64):
			sim_dets[i,:,k*64:(k+1)*64] -= np.random.choice(col_meds)
			
	sim_dets = sim_dets * np.sqrt(NDIT)
	return sim_dets


def add_detectors(cube, dets):
	''' Adds detectors to a datacube
	Inputs:
		cube: Input datacube (RA, DEC, lambda)
		dets: List of 8 detector arrays

	Outputs:
		new_datacube: Cube with detector systematics added
	'''
	sim_dets = np.copy(dets)
	cube_shape = cube.shape
	
	# allow for datacube to be smaller than the detectors
	full_datacube = np.zeros((3700,152,204))
	full_datacube[:cube_shape[0],:cube_shape[1],:cube_shape[2]] += cube

	# slice datacube up into 8 octants
	slice1 = full_datacube[:,:38,:102]
	slice2 = full_datacube[:,:38,102:]
	slice3 = full_datacube[:,38:76,:102]
	slice4 = full_datacube[:,38:76,102:]
	slice5 = full_datacube[:,76:114,:102]
	slice6 = full_datacube[:,76:114,102:]
	slice7 = full_datacube[:,114:,:102]
	slice8 = full_datacube[:,114:,102:]
	slices = [slice1, slice2, slice3, slice4, slice5, slice6, slice7, slice8]

		# add detector read noise onto datacube
	for i in range(len(sim_dets)):
		for j in range(slices[0].shape[1]):
			slitlet = slices[i][:,j,:]
			if i % 2 == 0:
				hoz_pos = np.int(np.round(32.5+(j*102)+(4.73*j)))
				if j % 2 == 0:
					vert_pos = 130 + (2 * j) - 67
				else:
					vert_pos = 130 + (2 * j) + 67
			else:
				hoz_pos = np.int(np.round(12.5+(j*102)+(4.73*j)))
				if j % 2 == 0:
					vert_pos = 204 - (2 * j) - 67
				else:
					vert_pos = 204 - (2 * j) + 67

			sim_dets[i][vert_pos:vert_pos+3700, hoz_pos:hoz_pos+102]+=slitlet

		# Convert detectors back into datacube		  
	new_datacube_slices=[]
	for i in range(len(sim_dets)):
		new_datacube_slice = np.zeros(slices[i].shape)
		for j in range(new_datacube_slice.shape[1]):
			if i % 2 == 0:
				hoz_pos = np.int(np.round(32.5+(j*102)+(4.73*j)))
				if j % 2 == 0:
					vert_pos = 130 + (2 * j) - 67
				else:
					vert_pos = 130 + (2 * j) + 67
			else:
				hoz_pos = np.int(np.round(12.5+(j*102)+(4.73*j)))
				if j % 2 == 0:
					vert_pos = 204 - (2 * j) - 67
				else:
					vert_pos = 204 - (2 * j) + 67
			new_datacube_slice[:,j,:] = sim_dets[i][vert_pos:vert_pos+3700, hoz_pos:hoz_pos+102]
		new_datacube_slices.append(new_datacube_slice)

	new_datacube = np.zeros((3700,152,204))

	new_datacube[:,:38,:102] = new_datacube_slices[0]
	new_datacube[:,:38,102:] = new_datacube_slices[1]
	new_datacube[:,38:76,:102] = new_datacube_slices[2]
	new_datacube[:,38:76,102:] = new_datacube_slices[3]
	new_datacube[:,76:114,:102] = new_datacube_slices[4]
	new_datacube[:,76:114,102:] = new_datacube_slices[5]
	new_datacube[:,114:,:102] = new_datacube_slices[6]
	new_datacube[:,114:,102:] = new_datacube_slices[7]

	# ensure new cube is same size as input cube
	new_datacube = new_datacube[:cube_shape[0],:cube_shape[1],:cube_shape[2]]
	return new_datacube
