# coding=utf-8
'''
Calculates PSF, jitter and telescope background and transmission
'''
import os
import sys
import logging

import multiprocessing as mp

import numpy as np
import scipy.constants as sp
from scipy.interpolate import interp1d
from scipy.signal import fftconvolve

from src.config import *

from src.modules.create_psf import create_psf, define_psf
from src.modules.misc_utils import path_setup
from src.modules.blackbody import *

import matplotlib.pylab as plt

tppath = path_setup('../../' + config_data["data_dir"] + 'throughput/')


def telescope_background_emission(wavels, T, emissivity, DIT, debug_plots, output_file):
	'''Function that generates a telescope background curve
	using mirror reflectivities, mirror areas and the
	Plank BB function as a function of site temperature.
	Telescope modelled as a graybody with eps*BB(T), where
	eps = emissivity of telescope. Emissivity is computed with
	mirror reflectivites for 6 mirrors of telescope.
	
	Inputs:
		wavels: array of wavelengths for datacube
		T: site temperature [K]
		emissivity: Telescope emissivity modelled as 1 - reflectivity
		dit: exposure time [s]. This determins how the sky emission
		
	Outputs:
		tele_bg_curve: array of total telescope background
			[units of photons/m^2/um/arcsec^2]
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

	tele_bg_spec = np.multiply(emissivity, cube_bb_spec)

	tele_bg_spec_ph = tele_bg_spec/(sp.h*sp.c/(wavels*1.E-6))*DIT

	if debug_plots:
		plt.clf()
		plt.plot(wavels, tele_bg_spec_ph, label="Blackbody T = {:.1f} K".format(T))
		plt.legend()
		plt.xlabel(r"wavelength [$\mu$m]")
		plt.ylabel(r"telescope emission [photons/m$^2$/$\mu$m/arcsec$^2$]")
		plt.savefig(output_file + "_tel_em.pdf")
		np.savetxt(output_file + "_tel_em.txt", np.c_[wavels, tele_bg_spec_ph])
	
	return tele_bg_spec_ph



#Telescope throughput curve generated just using wavelength array.
def telescope_transmission_curve(wavels, debug_plots, output_file):
	'''Function that reads a full telescope throughput curve.

	Inputs:
		wavels: array of wavelengths for datacube

	Outputs:
		cube_tele_trans: array of telescope throughput
			for each wavelength value in wavels
	'''

	#Load telescope reflectivity file
	tele_r = np.genfromtxt(os.path.join(tppath, 'ELT_mirror_reflectivity.txt'), delimiter=',')

	#Interpolate as a function of wavelength
	tele_trans_interp = interp1d(tele_r[:, 0], tele_r[:, 1],
				kind='linear', bounds_error=False, fill_value=0.)
	#Obtain values for datacube wavelength array
	cube_tele_trans = tele_trans_interp(wavels)
	
	
	if debug_plots:
		plt.clf()
		plt.plot(wavels, cube_tele_trans)
		plt.xlabel(r"wavelength [$\mu$m]")
		plt.ylabel(r"telescope transmission ")
		plt.savefig(output_file + "_tel_tr.pdf")
		np.savetxt(output_file + "_tel_tr.txt", np.c_[wavels, cube_tele_trans])



	return cube_tele_trans


def process_lambda(params, lamb, image, px, py, pback):
	padding_plane = np.zeros((px, py))
	#return i, padding_plane
	psf = create_psf(lamb)
	# np.sum(psf) is < 1, so to recover the same back emission after the convolution we need to 
	# subctract before and then add the back emission. We assume that the back emission is uniform
	# within the field
	padding_plane[psf.shape[0]:psf.shape[0]+image.shape[0], psf.shape[1]:psf.shape[1]+image.shape[1]] = image - pback
	conv_image = fftconvolve(padding_plane, psf) + pback
	return params, conv_image
	
	
counter = None
result_cube = None
bar_str = None
llambs = None

def save_result(results):
	global counter, result_cube, bar_str, llambs
	counter = counter + 1
	sys.stdout.write(bar_str.format(int(100.*counter/llambs), counter, llambs) + "\r")
	sys.stdout.flush()

	(i, x0, x1, y0, y1), conv_image = results
	result_cube[i, :, :] = conv_image[y0:y1, x0:x1]


def sim_telescope(cube, back_emission, ext_lambs, cube_lamb_mask, DIT, jitter, air_mass, seeing, spax, site_temp, aoMode, ncpus, debug_plots=False, output_file=""):
	''' Simulates telescope effects
	Inputs:
		cube: Input datacube (RA, DEC, lambda)
		back_emission: Input background emission outside of the FoV
		ext_lambs: extended lambda array [um]
		cube_lamb_mask: mask array to get the lambs of the cube
		DIT: Exposure time [s]
		jitter: Residual telescope jitter [mas]
		air_mass: Air mass of the observation
		seeing: Atmospheric seeing FWHM [arcsec]
		spax: spatial pixel (spaxel) scale [mas]
		site_temp: Telescope temperature [K]
		aoMode: AO mode: LTAO/SCAO/NOAO/AIRY
		ncpus: no. of CPUs to use
		debug_plots: Produce debug plots
		output_file: File name for debug plots
		
	Outputs:
		cube: Cube including telescope background and emission and convolved with PSF
		back_emission: back_emission including telescope
		PSF: PSF of the first lambda
	'''
	
	# Get telescope reflectivity
	logging.info("Calculating telescope reflectivity")
	telescope_reflectivity = telescope_transmission_curve(ext_lambs, debug_plots, output_file)
	back_emission = back_emission*telescope_reflectivity
	
	# Get telescope background
	logging.info("Calculating telescope background")
	telescope_background = telescope_background_emission(ext_lambs, site_temp, 1. - telescope_reflectivity, DIT, debug_plots, output_file)
	back_emission = back_emission + telescope_background
	
	# Add telescope emission/transmission to the input cube
	tel_reflectivity_cube = telescope_reflectivity[cube_lamb_mask]
	tel_reflectivity_cube.shape = (np.sum(cube_lamb_mask), 1, 1)
	cube *= tel_reflectivity_cube

	tel_background_cube = telescope_background[cube_lamb_mask]
	tel_background_cube.shape = (np.sum(cube_lamb_mask), 1, 1)
	cube += tel_background_cube
	
	
	# PSF + Jitter + Instrument PSF
	logging.info("Define PSF")
	
	FWHM_instrument = (config_data["dynamic_instrument_psf"]**2 + config_data["static_instrument_psf"][spax]**2)**0.5
	sigma_instrument = FWHM_instrument/2.35482
	sigma_combined = (jitter**2 + sigma_instrument**2)**0.5
	
	logging.info("Residual telescope jitter = {:.2f} mas".format(jitter))
	logging.info("Instrument PSF = {:.2f} mas".format(sigma_instrument))
	logging.info("-> Combined σ = {:.2f} mas".format(sigma_combined))
	
	psf_size = config_data["spaxel_scale"][spax].psfsize
	
	define_psf(air_mass, seeing, sigma_combined, config_data["telescope"]["diameter"], psf_size, config_data["spaxel_scale"][spax].psfscale, aoMode)
	lambs = ext_lambs[cube_lamb_mask]
	
	# padding with back_emission
	padding_x = cube.shape[2] + 2*psf_size
	padding_y = cube.shape[1] + 2*psf_size
	
	padding_back = back_emission[cube_lamb_mask]
	padding_back.shape = (len(lambs), 1, 1)
	
	conv_side_x = padding_x + psf_size - 1
	conv_side_y = padding_y + psf_size - 1

	# extract convolved image
	x0 = int((conv_side_x - cube.shape[2])/2.)+1
	x1 = x0 + cube.shape[2]
	y0 = int((conv_side_y - cube.shape[1])/2.)+1
	y1 = y0 + cube.shape[1]
	
	logging.info("Convolve with PSF")
	
	global counter, result_cube, bar_str, llambs
	
	bar_str = "[ {:2d}% {:" + str(len(str(len(lambs)))) + "d}/{:" + str(len(str(len(lambs)))) + "d} ]"
	
	logging.info(("Using {:d} CPU" + ("s" if ncpus > 1 else "")).format(ncpus))
	if ncpus > 1:

		pool = mp.Pool(processes=ncpus)
		
		counter = 0
		result_cube = np.zeros_like(cube)
		llambs = len(lambs)
		
		for i in range(len(lambs)):
			pool.apply_async(process_lambda, args=((i, x0, x1, y0, y1), lambs[i], cube[i, :, :], padding_x, padding_y, padding_back[i]), callback=save_result)
		
		pool.close()
		pool.join()
		
		cube = result_cube

	else:
		for i in range(len(lambs)):
			sys.stdout.write(bar_str.format(int(100.*i/len(lambs)), i, len(lambs)) + "\r")
			sys.stdout.flush()
			
			_, conv_image = process_lambda(i, lambs[i], cube[i, :, :], padding_x, padding_y, padding_back[i])
		
			cube[i, :, :] = conv_image[y0:y1, x0:x1]
			#break
	
	central_lambda = (lambs[0] + lambs[-1])*0.5
	return cube, back_emission, create_psf(central_lambda), central_lambda


