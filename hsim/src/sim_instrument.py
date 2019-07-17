'''
Calculates LSF, instrument background and transmission
'''
import logging

import numpy as np
import scipy.constants as sp
from astropy.convolution import Gaussian1DKernel

from src.config import *
from src.modules.em_model import *


def sim_instrument(cube, back_emission, transmission, ext_lambs, cube_lamb_mask, DIT, grating, site_temp, input_spec_res, aoMode, debug_plots=False, output_file=""):
	''' Simulates instrument effects
	Inputs:
		cube: Input datacube (RA, DEC, lambda)
		back_emission: Input background emission
		transmission: Input transmission
		ext_lambs: extended lambda array [um]
		cube_lamb_mask: mask array to get the lambs of the cube
		DIT: Exposure time [s]
		grating: Spectral grating
		site_temp: Telescope temperature [K]
		input_spec_res: Spectral resolution of the input cube [micron]
		aoMode: LTAO/SCAO/NOAO
		debug_plots: Produce debug plots
		output_file: File name for debug plots

	Outputs:
		cube: cube including instrument effects
		back_emission: back_emission including telescope
		LSF_size: width of the LSF [A]
	'''
	
	# Get instrument transmission
	logging.info("Calculating HARMONI transmission")
	# LTAO dichroic if present
	if aoMode in ["LTAO", "SCAO"]:
		AOd_tr = load_transmission_curve(ext_lambs, "LTAO_dichroic.txt", debug_plots, [output_file, "ins_AOd"], "AO dichroic transmission")
	else:
		logging.warning("AO Dichroic throughput file not found for " + aoMode)
		AOd_tr = 1.
	
	# FPRS
	FPRS_tr = load_transmission_curve(ext_lambs, "FPRS.txt", debug_plots, [output_file, "ins_FPRS"], "FPRS transmission")
	
	
	# instrument after FPRS
	# Cryostat, Pre-optics, IFU, Spectrograph and Grating
	instrument_tr = load_transmission_curve(ext_lambs, grating + '_grating.txt', debug_plots, [output_file, "ins"], "instrument transmission")
	

	back_emission = back_emission*AOd_tr*FPRS_tr*instrument_tr
	transmission = transmission*AOd_tr*FPRS_tr*instrument_tr
	
	# Get instrument background
	logging.info("Calculating HARMONI background")

	# AO dichroic if present
	if aoMode in ["LTAO", "SCAO"]:
		AOd_background = get_background_emission(ext_lambs, site_temp - config_data['HARMONI_AO_diff_temp'][aoMode], 1. - AOd_tr, DIT, debug_plots, [output_file, "ins", "AOd"], "instrument emission [photons/m$^2$/$\mu$m/arcsec$^2$]")
	else:
		AOd_background = 0.

	# FPRS
	FPRS_background = get_background_emission(ext_lambs, site_temp - config_data['HARMONI_FPRS_diff_temp'], 1. - FPRS_tr, DIT, debug_plots, [output_file, "ins", "FPRS"], "instrument emission [photons/m$^2$/$\mu$m/arcsec$^2$]")
		
	# instrument after FPRS
	instrument_background = get_background_emission(ext_lambs, config_data['HARMONI_temp'], 1. - instrument_tr, DIT, debug_plots, [output_file, "ins"], "instrument emission [photons/m$^2$/$\mu$m/arcsec$^2$]")
	
	back_emission = back_emission + AOd_background + FPRS_background + instrument_background
	
	# Add instrument emission/transmission to the input cube
	
	instrument_tr_cube = (AOd_tr*FPRS_tr*instrument_tr)[cube_lamb_mask]
	instrument_tr_cube.shape = (np.sum(cube_lamb_mask), 1, 1)
	cube *= instrument_tr_cube

	total_instrument_background = AOd_background + FPRS_background + instrument_background
	instrument_background_cube = total_instrument_background[cube_lamb_mask]
	instrument_background_cube.shape = (np.sum(cube_lamb_mask), 1, 1)
	cube += instrument_background_cube
	
	
	logging.info("Convolve with LSF")
	# Assume Gaussian LSF
	bandws = config_data['gratings'][grating]
	new_res = (bandws.lmin + bandws.lmax)/(2.*bandws.R) # micron
	new_res_pix = (new_res**2 - input_spec_res**2)**0.5/(ext_lambs[1] - ext_lambs[0])
	logging.info("Output resolution: {:.3f} A".format(new_res*10000.))
	logging.info("Input resolution: {:.3f} A".format(input_spec_res*10000.))
	logging.info("LSF FWHM = {:.3f} A".format((new_res**2 - input_spec_res**2)**0.5*10000.))
	
	
	npix_LSF = 0
	if new_res_pix > 0.:
		sigma_LSF_pix = new_res_pix/2.35482
	
		npix_LSF = int(sigma_LSF_pix*config_data['LSF_kernel_size'])
		kernel_LSF = Gaussian1DKernel(stddev=sigma_LSF_pix, x_size=npix_LSF)
		z, y, x = cube.shape
		
		for py in range(y):
			for px in range(x):
				spectrum = np.copy(back_emission)
				spectrum[cube_lamb_mask] = cube[:, py, px]
				
				cube[:, py, px] = np.convolve(spectrum, kernel_LSF, mode="same")[cube_lamb_mask]
		
		
		back_emission = np.convolve(back_emission, kernel_LSF, mode="same")
		transmission = np.convolve(transmission, kernel_LSF, mode="same")

	LSF_size = npix_LSF*(ext_lambs[1] - ext_lambs[0])*10000. # A
	logging.info("Total LSF width for the convolution: {:.3f} A".format(LSF_size))
	
	return cube, back_emission, transmission, LSF_size
