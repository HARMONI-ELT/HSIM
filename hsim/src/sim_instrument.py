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


def HARMONI_background(wavels, Temp, DIT, debug_plots, output_file):
	'''Function that generates a HARMONI background curve
	using mirror reflectivities, mirror areas and the
	Plank BB function as a function of instrument temperature.
	HARMONI modelled as a blackbody BB(T) with T = 130 K.
	
	Inputs:
		wavels: array of wavelengths for datacube
		Temp: instrument temperature
		DIT: Exposure time [s]
		
	Outputs:
		inst_bg_spec_ph: array of total instrument background
			[units of photons/s/m^2/um/arcsec^2]
			for each wavelength value in wavels
	'''

	#Blackbody function
	cube_bb_spec = blackbody(wavels, Temp)

	inst_bg_spec = cube_bb_spec
	inst_bg_spec_ph = inst_bg_spec/(sp.h*sp.c/(wavels*1.E-6))*DIT
	
	if debug_plots:
		plt.clf()
		plt.plot(wavels, inst_bg_spec_ph, label="Blackbody T = {:.1f} K".format(Temp))
		plt.legend()
		plt.xlabel(r"wavelength [$\mu$m]")
		plt.ylabel(r"instrument emission [photons/m$^2$/$\mu$m/arcsec$^2$]")
		plt.savefig(output_file + "_ins_em.pdf")
		np.savetxt(output_file + "_ins_em.txt", np.c_[wavels, inst_bg_spec_ph])
	
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
	inst_r = np.genfromtxt(os.path.join(tppath, grating+'_grating.txt'), delimiter=',')

	#Interpolate as a function of wavelength
	inst_trans_interp = interp1d(inst_r[:,0], inst_r[:,1],
					kind='linear', bounds_error=False, fill_value=0.)
	
	#Obtain values for datacube wavelength array
	cube_inst_trans = inst_trans_interp(wavels)

	if debug_plots:
		plt.clf()
		plt.plot(wavels, cube_inst_trans)
		plt.xlabel(r"wavelength [$\mu$m]")
		plt.ylabel(r"instrument transmission ")
		plt.savefig(output_file + "_ins_tr.pdf")
		np.savetxt(output_file + "_ins_tr.txt", np.c_[wavels, cube_inst_trans])


	return cube_inst_trans



def sim_instrument(cube, back_emission, ext_lambs, cube_lamb_mask, DIT, grating, debug_plots=False, output_file=""):
	''' Simulates instrument effects
	Inputs:
		cube: Input datacube (RA, DEC, lambda)
		back_emission: Input background emission outside of the FoV
		ext_lambs: extended lambda array [um]
		cube_lamb_mask: mask array to get the lambs of the cube
		DIT: Exposure time [s]
		grating: Spectral grating
		debug_plots: Produce debug plots
		output_file: File name for debug plots

	Outputs:
		cube: Cube including ...
		back_emission: back_emission including telescope
	'''
	
	# Get instrument transmission
	print "Calculating HARMONI transmission"
	instrument_tr = HARMONI_transmission_curve(ext_lambs, grating, debug_plots, output_file)
	back_emission = back_emission*instrument_tr
	
	# Get instrument background
	print "Calculating HARMONI background"
	instrument_background = HARMONI_background(ext_lambs, config_data['HARMONI_temp'], DIT, debug_plots, output_file)
	back_emission = back_emission + instrument_background
	
	# Add instrument emission/transmission to the input cube
	
	instrument_tr_cube = instrument_tr[cube_lamb_mask]
	instrument_tr_cube.shape = (np.sum(cube_lamb_mask),1,1)
	cube = np.multiply(cube, instrument_tr_cube)

	instrument_background_cube = instrument_background[cube_lamb_mask]
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
	

	