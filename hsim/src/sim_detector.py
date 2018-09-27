'''
Calculates detector QE, dark current, read noise
'''
import os
import logging

import numpy as np
import scipy.constants as sp
from scipy.interpolate import interp1d
from astropy.convolution import Gaussian1DKernel

from config import *

from modules.misc_utils import path_setup
from modules.blackbody import *

import matplotlib.pylab as plt

tppath = path_setup('../../' + config_data["data_dir"] + 'throughput/')

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
		
	det_qe = np.genfromtxt(os.path.join(tppath, detector_QE_file), delimiter=',')


	#Interpolate as a function of wavelength
	det_qe_interp = interp1d(det_qe[:,0], det_qe[:,1],
				kind='linear', bounds_error=False, fill_value=0.)
	#Obtain values for datacube wavelength array
	cube_det_qe = det_qe_interp(wavels)/100.


	if debug_plots:
		plt.clf()
		plt.plot(wavels, cube_det_qe)
		plt.xlabel(r"wavelength [$\mu$m]")
		plt.ylabel(r"detector QE ")
		plt.savefig(output_file + "_det_qe.pdf")
		np.savetxt(output_file + "_det_qe.txt", np.c_[wavels, cube_det_qe])


	return cube_det_qe


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
	
	logging.info("Applying detector crosstalk")
	
	scaled_cube = cube*(1. - crosstalk*4)
	
	# crosstalk in the spatial direction
	spatial_crosstalk = crosstalk*(np.roll(cube, 1, axis=2) + np.roll(cube, -1, axis=2))
	
	# crosstalk in the spectral direction
	spectral_crosstalk = crosstalk*(np.roll(cube, 1, axis=0) + np.roll(cube, -1, axis=0))
	
	return scaled_cube + spatial_crosstalk + spectral_crosstalk

 

def sim_detector(cube, back_emission, lambs, grating, DIT, debug_plots=False, output_file=""):
	''' Simulates detector effects
	Inputs:
		cube: Input datacube (RA, DEC, lambda)
		back_emission: Input background emission
		lambs: lambda array [um]
		grating: Spectral grating
		debug_plots: Produce debug plots
		output_file: File name for debug plots

	Outputs:
		cube: Cube including QE
		read_noise: read noise for the grating and DIT [e/pix]
		dark*DIT: dark current  [e/pix]
	'''
	
	# Get QE curve
	logging.info("Calculating detector QE")
	qe_curve = detector_QE_curve(lambs, grating, debug_plots, output_file)
	back_emission = np.multiply(back_emission, qe_curve)
	
	qe_curve.shape = (len(lambs),1,1)
	cube = np.multiply(cube, qe_curve)

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

	
	return cube, back_emission, read_noise, dark*DIT
	

	