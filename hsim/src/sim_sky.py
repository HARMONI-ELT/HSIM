'''
Calculates the sky transmission and emission at the observed lambda
and computes the ADR if requested
'''
import os
import logging

import numpy as np
from astropy.convolution import Gaussian1DKernel

from src.config import *
from src.modules.misc_utils import path_setup
from src.modules.rebin import *
from src.modules.adr import apply_adr

import matplotlib.pyplot as plt

bgpath = path_setup('../../' + config_data["data_dir"] + 'sky/')

def convolve_1d_spectrum(input_lambda, input_flux, output_spec_res):
	'''Function that convolves a sky spectrum with a Gaussian to 
	match the input cube spectral resolution.
	
	Inputs:
		input_lambda: input sky spectrum lambda
		input_flux: input sky spectrum flux
		output_spec_res: Spectral resolution of the convolved sky spectrum output 
		
	Outputs:
		convolved sky spectrum
	'''
	
	sky_resolution = np.abs(input_lambda[1] - input_lambda[0])
	
	if output_spec_res > sky_resolution:
		logging.info("Convolve sky to input cube spectral resolution")
		
		new_res_pix = (output_spec_res**2 - sky_resolution**2)**0.5/np.abs(input_lambda[1] - input_lambda[0])
		sigma_LSF_pix = new_res_pix/2.35482
		if sigma_LSF_pix < 1.: # avoid convolution with a kernel narrower than 1 pixel
			return input_flux
		npix_LSF = int(sigma_LSF_pix*config_data['LSF_kernel_size'])
		kernel_LSF = Gaussian1DKernel(stddev=sigma_LSF_pix, x_size=npix_LSF)
		return np.convolve(input_flux, kernel_LSF, mode="same")
		
	else:
		logging.warning("Sky spectral resolution (R = {:.0f}) is lower than the input cube spectral resolution (R = {:.0f})".format(np.median(input_lambda)/sky_resolution, np.median(input_lambda)/output_spec_res if output_spec_res != 0 else np.inf))
		return input_flux
	
	
def sky_background(input_parameters, lambs, air_mass, dit, input_spec_res, debug_plots, output_file):
	'''Function that generates a sky background curve combining
	sky continuum, sky thermal emission and sky emission lines.
	
	Inputs:
		lambs: array of wavelengths for datacube
		air_mass: Air mass of the observation
		dit: exposure time [s].
		input_spec_res: Spectral resolution of the input cube [micron]
		
	Outputs:
		sky_radiance: array of total sky background for DIT
			[units of photons/m^2/um/arcsec^2]
	'''
	inbuilt_airmasses = [1.1, 1.3, 1.5, 2.0]
	if air_mass not in inbuilt_airmasses:
		raise HSIMError('Error: ' + str(air_mass) + ' is not a valid air_mass. Valid options are: ' + ",".join([str(_) for _ in inbuilt_airmasses]))


	#determine the closest data to the airmass value given and find it's location in the data file
	closest_X = min(inbuilt_airmasses, key=lambda x:abs(x - air_mass))
	data_index = inbuilt_airmasses.index(closest_X) + 1

	#load sky transmission & extinction files, then reduce to the columns required
	sky_em_all_X = np.genfromtxt(os.path.join(bgpath, 'radiance.txt.gz'), delimiter=',')
	
	sky_em_lambda = sky_em_all_X[:, 0]
	sky_em_flux = sky_em_all_X[:, data_index]
	
	# estimate scattered background
	bandws = config_data['gratings'][input_parameters["grating"]]
	mask_range_scatter = (sky_em_lambda > bandws.lmin)*(sky_em_lambda < bandws.lmax)
	mean_sky = np.mean(sky_em_flux[mask_range_scatter])
	if input_parameters["mci"]:
		scattered_sky = mean_sky*0.2
	else:
		scattered_sky = mean_sky*input_parameters["scattered_sky"]/100.
	
	logging.info("Extra scatter sky = " + str(scattered_sky) + " photons/m2/um/arcsec2")
	

	mask_range_output = (sky_em_lambda > lambs[0]*0.99)*(sky_em_lambda < lambs[-1]*1.01)
	sky_em_lambda = sky_em_lambda[mask_range_output]
	sky_em_flux = sky_em_flux[mask_range_output] + scattered_sky
	
	# Match input cube spectral resolution
	sky_em_flux = convolve_1d_spectrum(sky_em_lambda, sky_em_flux, input_spec_res)

	# rebin sky emission
	sky_radiance = dit*rebin1d(lambs, sky_em_lambda, sky_em_flux)

	if debug_plots:
		plt.clf()
		mask_plot = (sky_em_lambda > lambs[0])*(sky_em_lambda < lambs[-1])
		plt.plot(sky_em_lambda[mask_plot], dit*sky_em_flux[mask_plot], label="Skycalc 0.15A")
		plt.plot(lambs, sky_radiance, label="rebin")
		plt.legend()
		plt.xlabel(r"wavelength [$\mu$m]")
		plt.ylabel(r"sky emission [photons/m$^2$/$\mu$m/arcsec$^2$]")
		plt.savefig(output_file + "_sky_em.pdf")
		np.savetxt(output_file + "_sky_em.txt", np.c_[lambs, sky_radiance])


	
	return sky_radiance


def moon_background(lambs, moon, dit, input_spec_res, debug_plots, output_file):
	'''Function that generates a moon background curve
	
	Inputs:
		lambs: array of wavelengths for datacube
		moon: Fractional moon illumination
		dit: exposure time [s]. This determins how the sky emission
		input_spec_res: Spectral resolution of the input cube [micron]
		
	Outputs:
		sky_radiance: array of total sky background for DIT
			[units of photons/m^2/um/arcsec^2]
	'''
	
	if moon not in [0., 0.5, 1.0]:
		raise HSIMError('Error: ' + str(moon) + ' is not a valid Moon illumination. Valid options are: 0, 0.5, 1.0')
	
	if moon > 0.:
		inbuilt_moon = [0.5, 1.0]

		#determine the closest data to the airmass value given and find it's location in the data file
		closest_X = min(inbuilt_moon, key=lambda x: abs(x - moon))
		data_index = inbuilt_moon.index(closest_X) + 1

		#load sky transmission & extinction files, then reduce to the columns required
		moon_em_all_X = np.genfromtxt(os.path.join(bgpath, 'moon.txt.gz'), delimiter=',')
		
		moon_em_lambda = moon_em_all_X[:, 0]
		moon_em_flux = moon_em_all_X[:, data_index]

		mask_range_output = (moon_em_lambda > lambs[0]*0.99)*(moon_em_lambda < lambs[-1]*1.01)
		moon_em_lambda = moon_em_lambda[mask_range_output]
		moon_em_flux = moon_em_flux[mask_range_output]

		# Match input cube spectral resolution
		moon_em_flux = convolve_1d_spectrum(moon_em_lambda, moon_em_flux, input_spec_res)

		# rebin moon emission
		moon_radiance = dit*rebin1d(lambs, moon_em_lambda, moon_em_flux)
	else:
		moon_em_lambda = lambs
		moon_em_flux = lambs*0.
		moon_radiance = moon_em_flux
	
	if debug_plots:
		plt.clf()
		mask_plot = (moon_em_lambda > lambs[0])*(moon_em_lambda < lambs[-1])
		plt.plot(moon_em_lambda[mask_plot], dit*moon_em_flux[mask_plot], label="Skycalc 0.15A")
		plt.plot(lambs, moon_radiance, label="rebin")
		plt.legend()
		plt.xlabel(r"wavelength [$\mu$m]")
		plt.ylabel(r"moon emission [photons/m$^2$/$\mu$m/arcsec$^2$]")
		plt.savefig(output_file + "_moon_em.pdf")
		np.savetxt(output_file + "_moon_em.txt", np.c_[lambs, moon_radiance])

	return moon_radiance


#Sky throughput curve generated just using wavelength array.
def sky_transmission(lambs, air_mass, input_spec_res, debug_plots, output_file):
	'''Function that generates a full throughput curve combining
	sky transmission & sky extinction.

	Inputs:
		lambs: array of wavelengths for datacube
		air_mass: Air mass of the observation
		input_spec_res: Spectral resolution of the input cube [micron]

	Outputs:
		cube_total_sky_trans: array of total throughput
			for each wavelength value in lambs
	'''
	#convert from zenith angle to airmass
	inbuilt_airmasses = [1.1, 1.3, 1.5, 2.0]
	if air_mass not in inbuilt_airmasses:
		raise HSIMError('Error: ' + str(air_mass) + ' is not a valid air_mass. Valid options are: ' + ",".join([str(_) for _ in inbuilt_airmasses]))

	#determine the closest data to the airmass value given and find it's location in the data file
	closest_X = min(inbuilt_airmasses, key=lambda x: abs(x - air_mass))
	data_index = inbuilt_airmasses.index(closest_X) + 1

	#load sky transmission & extinction files, then reduce to the columns required
	sky_trans_all_X = np.genfromtxt(os.path.join(bgpath, 'transmission.txt.gz'), delimiter=',')
	
	sky_tr_lambda = sky_trans_all_X[:, 0]
	sky_tr = sky_trans_all_X[:, data_index]

	mask_range_output = (sky_tr_lambda > lambs[0]*0.99)*(sky_tr_lambda < lambs[-1]*1.01)
	sky_tr_lambda = sky_tr_lambda[mask_range_output]
	sky_tr = sky_tr[mask_range_output]

	# Match input cube spectral resolution
	sky_tr = convolve_1d_spectrum(sky_tr_lambda, sky_tr, input_spec_res)

	final_tr = rebin1d(lambs, sky_tr_lambda, sky_tr)

	if debug_plots:
		plt.clf()
		mask_plot = (sky_tr_lambda > lambs[0])*(sky_tr_lambda < lambs[-1])
		plt.plot(sky_tr_lambda[mask_plot], sky_tr[mask_plot], label="Skycalc 0.15A")
		plt.plot(lambs, final_tr, label="rebin")
		plt.legend()
		plt.xlabel(r"wavelength [$\mu$m]")
		plt.ylabel(r"sky transmission")
		plt.savefig(output_file + "_sky_tr.pdf")
		np.savetxt(output_file + "_sky_tr.txt", np.c_[lambs, final_tr])

	return final_tr
	



def sim_sky(input_parameters, cube, back_emission, transmission, header, ext_lambs, cube_lamb_mask, input_spec_res, debug_plots=False, output_file=""):
	''' Simulates sky effects
	Inputs:
		input_parameters: input dictionary
			exposure_time: Exposure time [s]
			air_mass: Air mass of the observation
			moon_illumination: Fractional moon illumination
			telescope_temp: Telescope temperature [K]
			adr: Boolean - turn ADR on or off
		
		
		cube: Input datacube (RA, DEC, lambda)
		back_emission: Input background emission
		transmission: Input transmission
		header: FITS header
		ext_lambs: extended lambda array [um]
		cube_lamb_mask: mask array to get the lambs of the cube
		input_spec_res: Spectral resolution of the input cube [micron]
		debug_plots: Produce debug plots
		output_file: File name for debug plots
	Outputs:
		cube: Cube including sky emission, transmission and ADR
		back_emission: back_emission including sky
	'''
	air_mass = input_parameters["air_mass"]
	exposure_time = input_parameters["exposure_time"]
	moon = input_parameters["moon_illumination"]
	telescope_temp = input_parameters["telescope_temp"]
	adr_switch = input_parameters["adr"]
	
	# Get sky transmission
	logging.info("Calculating sky transmission")
	sky_trans = sky_transmission(ext_lambs, air_mass, input_spec_res, debug_plots, output_file)
	transmission *= sky_trans
	
	# Get sky emission (lines + continuum)
	logging.info("Calculating sky emission")
	sky_emission = sky_background(input_parameters, ext_lambs, air_mass, exposure_time, input_spec_res, debug_plots, output_file)
	
	# Get moon emission
	logging.info("Calculating Moon emission")
	moon_emission = moon_background(ext_lambs, moon, exposure_time, input_spec_res, debug_plots, output_file)
	back_emission = back_emission + sky_emission + moon_emission
	

	# Add sky emission/transmission to the input cube
	sky_trans_cube = sky_trans[cube_lamb_mask]
	sky_trans_cube.shape = (np.sum(cube_lamb_mask), 1, 1)
	cube *= sky_trans_cube

	sky_emission_cube = sky_emission[cube_lamb_mask] + moon_emission[cube_lamb_mask]
	sky_emission_cube.shape = (np.sum(cube_lamb_mask), 1, 1)
	cube += sky_emission_cube
	
	# Add atmospheric differential refration
	if adr_switch == True:
		logging.info("Calculating ADR")
		lambs = ext_lambs[cube_lamb_mask]
		cube = apply_adr(cube, header, lambs, telescope_temp, air_mass, debug_plots=False, output_file=output_file)
		
		
	return (cube, back_emission, transmission)

