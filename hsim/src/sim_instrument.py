'''
Calculates LSF, instrument background and transmission
'''
import os

import numpy as np
import scipy.constants as sp
from scipy.interpolate import interp1d
from astropy.convolution import Gaussian1DKernel

from config import *

from modules.misc_utils import path_setup
from modules.blackbody import *

import matplotlib.pylab as plt

tppath = path_setup('../../' + config_data["data_dir"] + 'throughput/')


def HARMONI_background(wavels, Temp, emissivity, element, DIT, debug_plots, output_file):
	'''Function that generates a HARMONI background curve
	using the Plank BB function as a function of instrument temperature.
	
	Inputs:
		wavels: array of wavelengths for datacube
		Temp: instrument temperature
		element: string to identify the HARMONI element for which the emission is calculated
		emissivity: emissivity
		DIT: Exposure time [s]
		
	Outputs:
		inst_bg_spec_ph: array of total instrument background
			[units of photons/m^2/um/arcsec^2]
			for each wavelength value in wavels
	'''

	#Blackbody function
	cube_bb_spec = blackbody(wavels, Temp)

	inst_bg_spec = np.multiply(emissivity, cube_bb_spec)
	
	inst_bg_spec_ph = inst_bg_spec/(sp.h*sp.c/(wavels*1.E-6))*DIT
	
	if debug_plots:
		plt.clf()
		plt.plot(wavels, inst_bg_spec_ph, label="Blackbody T = {:.1f} K".format(Temp))
		plt.legend()
		plt.xlabel(r"wavelength [$\mu$m]")
		plt.ylabel(r"instrument emission [photons/m$^2$/$\mu$m/arcsec$^2$]")
		plt.savefig(output_file + "_ins_"  + (element + "_" if element != "" else "") + "em.pdf")
		np.savetxt(output_file + "_ins_" + (element + "_" if element != "" else "") + "em.txt", np.c_[wavels, inst_bg_spec_ph])
	
	return inst_bg_spec_ph


#Instrument throughput curve generated just using wavelength array.
def HARMONI_transmission_curve(wavels, grating, debug_plots, output_file):
	'''Function that generates a full HARMONI throughput curve.

	Inputs:
		wavels: array of wavelengths for datacube
		grating: grating choice to set grating throughput curve

	Outputs:
		cube_inst_trans: array of instrument throughput
			for each wavelength value in wavels
	'''

	#Load grating throughput file
	# Cryostat, Pre-optics, IFU, Spectrograph and Grating
	inst_r = np.genfromtxt(os.path.join(tppath, grating+'_grating.txt'), delimiter=',')

	#Interpolate as a function of wavelength
	inst_trans_interp = interp1d(inst_r[:,0], inst_r[:,1],
					kind='linear', bounds_error=False, fill_value=0.)
	
	#Obtain values for datacube wavelength array
	cube_inst_trans = inst_trans_interp(wavels)

	if debug_plots:
		plt.clf()
		plt.plot(wavels, cube_inst_trans, label="after FPRS")
		plt.xlabel(r"wavelength [$\mu$m]")
		plt.ylabel(r"instrument transmission ")
		plt.legend(loc=2)
		plt.savefig(output_file + "_ins_tr.pdf")
		np.savetxt(output_file + "_ins_tr.txt", np.c_[wavels, cube_inst_trans])

	
	return cube_inst_trans


def LTAO_dichroic_transmission_curve(wavels, grating, debug_plots, output_file):
	'''Function that generates a full HARMONI throughput curve.

	Inputs:
		wavels: array of wavelengths for datacube
		grating: grating choice to set grating throughput curve

	Outputs:
		cube_ltaod_trans: array of LTAO dichroic throughput
			for each wavelength value in wavels
	'''

	#Load LTAO dichroic throughput file
	inst_r = np.genfromtxt(os.path.join(tppath, "LTAO_dichroic.txt"), delimiter=',')

	#Interpolate as a function of wavelength
	ltaod_trans_interp = interp1d(inst_r[:,0], inst_r[:,1],
					kind='linear', bounds_error=False, fill_value=0.)
	
	#Obtain values for datacube wavelength array
	cube_ltaod_trans = ltaod_trans_interp(wavels)

	if debug_plots:
		plt.clf()
		plt.plot(wavels, cube_ltaod_trans)
		plt.xlabel(r"wavelength [$\mu$m]")
		plt.ylabel(r"instrument transmission ")
		plt.savefig(output_file + "_ins_LTAOd_tr.pdf")
		np.savetxt(output_file + "_ins_LTAOd_tr.txt", np.c_[wavels, cube_ltaod_trans])

	
	
	return cube_ltaod_trans


def FPRS_transmission_curve(wavels, grating, debug_plots, output_file):
	'''Function that generates a full HARMONI throughput curve.

	Inputs:
		wavels: array of wavelengths for datacube
		grating: grating choice to set grating throughput curve

	Outputs:
		cube_fprs_trans: array of FPRS throughput
			for each wavelength value in wavels
	'''

	#Load FPRS throughput file
	inst_r = np.genfromtxt(os.path.join(tppath, "FPRS.txt"), delimiter=',')

	#Interpolate as a function of wavelength
	fprs_trans_interp = interp1d(inst_r[:,0], inst_r[:,1],
					kind='linear', bounds_error=False, fill_value=0.)
	
	#Obtain values for datacube wavelength array
	cube_fprs_trans = fprs_trans_interp(wavels)
	
	if debug_plots:
		plt.clf()
		plt.plot(wavels, cube_fprs_trans)
		plt.xlabel(r"wavelength [$\mu$m]")
		plt.ylabel(r"instrument transmission ")
		plt.savefig(output_file + "_ins_FPRS_tr.pdf")
		np.savetxt(output_file + "_ins_FPRS_tr.txt", np.c_[wavels, cube_fprs_trans])

	
	return cube_fprs_trans



def sim_instrument(cube, back_emission, ext_lambs, cube_lamb_mask, DIT, grating, site_temp, debug_plots=False, output_file="", LTAO_dichroic=True):
	''' Simulates instrument effects
	Inputs:
		cube: Input datacube (RA, DEC, lambda)
		back_emission: Input background emission outside of the FoV
		ext_lambs: extended lambda array [um]
		cube_lamb_mask: mask array to get the lambs of the cube
		DIT: Exposure time [s]
		grating: Spectral grating
		site_temp: Telescope temperature [K]
		debug_plots: Produce debug plots
		output_file: File name for debug plots
		LTAO_dichroic: True if the LTAO dichroic is used

	Outputs:
		cube: Cube including ...
		back_emission: back_emission including telescope
	'''
	
	# Get instrument transmission
	print "Calculating HARMONI transmission"
	# LTAO dichroic if present
	if LTAO_dichroic:
		LTAOd_tr = LTAO_dichroic_transmission_curve(ext_lambs, grating, debug_plots, output_file)
	else:
		LTAOd_tr = 1.
	
	# FPRS
	FPRS_tr = FPRS_transmission_curve(ext_lambs, grating, debug_plots, output_file)
	
	# instrument after FPRS
	instrument_tr = HARMONI_transmission_curve(ext_lambs, grating, debug_plots, output_file)

	back_emission = back_emission*LTAOd_tr*FPRS_tr*instrument_tr
	
	
	
	# Get instrument background
	print "Calculating HARMONI background"

	# LTAO dichroic if present
	if LTAO_dichroic:
		LTAOd_background = HARMONI_background(ext_lambs, site_temp - config_data['HARMONI_LTAO_diff_temp'], 1. - LTAOd_tr, "LTAOd", DIT, debug_plots, output_file)
	else:
		LTAOd_background = 0.

	# FPRS
	FPRS_background = HARMONI_background(ext_lambs, site_temp - config_data['HARMONI_FPRS_diff_temp'], 1. - FPRS_tr, "FPRS", DIT, debug_plots, output_file)
	
	# instrument after FPRS
	instrument_background = HARMONI_background(ext_lambs, config_data['HARMONI_temp'], 1. - instrument_tr, "", DIT, debug_plots, output_file)
	
	back_emission = back_emission + LTAOd_background + FPRS_background + instrument_background
	
	# Add instrument emission/transmission to the input cube
	
	instrument_tr_cube = (LTAOd_tr*FPRS_tr*instrument_tr)[cube_lamb_mask]
	instrument_tr_cube.shape = (np.sum(cube_lamb_mask),1,1)
	cube = np.multiply(cube, instrument_tr_cube)

	total_instrument_background = LTAOd_background + FPRS_background + instrument_background
	instrument_background_cube = total_instrument_background[cube_lamb_mask]
	instrument_background_cube.shape = (np.sum(cube_lamb_mask),1,1)
	cube = np.add(cube, instrument_background_cube)
	
	
	print "Convolve with LSF"
	# Assume Gaussian LSF
	bandws = config_data['gratings'][grating]
	new_res = (bandws.lmin + bandws.lmax)/(2.*bandws.R)
	new_res_pix = new_res/(ext_lambs[1] - ext_lambs[0])
	
	sigma_LSF_pix = new_res_pix/2.35482
	
	kernel_LSF = Gaussian1DKernel(stddev = sigma_LSF_pix, x_size = int(sigma_LSF_pix*config_data['LSF_kernel_size']))
	z, y, x = cube.shape
	
	for py in range(y):
		for px in range(x):
			spectrum = np.copy(back_emission)
			spectrum[cube_lamb_mask] = cube[:,py,px]
			
			cube[:,py,px] = np.convolve(spectrum, kernel_LSF, mode="same")[cube_lamb_mask]
	
	
	back_emission = np.convolve(back_emission, kernel_LSF, mode="same")

	return cube, back_emission
	

	