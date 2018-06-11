'''
Calculates the sky transmission and emission at the observed lambda
and computes the ADR if requested
'''
import os
import numpy as np

from config import *
from modules.misc_utils import path_setup
from modules.rebin import *
from modules.adr import apply_adr

import matplotlib.pylab as plt

bgpath = path_setup('../../' + config_data["data_dir"] + 'sky/')


def sky_background(lambs, air_mass, dit, debug_plots, output_file):
	'''Function that generates a sky background curve combining
	sky continuum, sky thermal emission and sky emission lines.
	
	Inputs:
		lambs: array of wavelengths for datacube
		air_mass: Air mass of the observation

		dit: exposure time [s]. This determins how the sky emission
		line amplitudes vary through the exposure.
		
	Outputs:
		sky_radiance: array of total sky background for DIT
			[units of photons/m^2/um/arcsec^2]
	'''
	inbuilt_airmasses = [1.1, 1.3, 1.5, 2.0]
	if air_mass not in inbuilt_airmasses:
		raise ValueError('Error: ' + str(air_mass) + ' is not a valid air_mass. Valid options are: ' + ",".join([str(_) for _ in inbuilt_airmasses]))


	#determine the closest data to the airmass value given and find it's location in the data file
	closest_X = min(inbuilt_airmasses, key=lambda x:abs(x - air_mass))
	data_index = inbuilt_airmasses.index(closest_X) + 1

	#load sky transmission & extinction files, then reduce to the columns required
	sky_em_all_X = np.genfromtxt(os.path.join(bgpath, 'radiance.txt'), delimiter=',')
	
	sky_em_lambda = sky_em_all_X[:,0]
	sky_em_flux = sky_em_all_X[:,data_index]

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


def moon_background(lambs, moon, dit, debug_plots, output_file):
	'''Function that generates a moon background curve
	
	Inputs:
		lambs: array of wavelengths for datacube
		moon: Fractional moon illumination

		dit: exposure time [s]. This determins how the sky emission
		
	Outputs:
		sky_radiance: array of total sky background for DIT
			[units of photons/m^2/um/arcsec^2]
	'''
	
	if moon not in [0., 0.5, 1.0]:
		raise ValueError('Error: ' + str(moon) + ' is not a valid Moon illumination. Valid options are: 0, 0.5, 1.0')
	
	if moon > 0.:
		inbuilt_moon = [0.5, 1.0]

		#determine the closest data to the airmass value given and find it's location in the data file
		closest_X = min(inbuilt_moon, key=lambda x:abs(x - moon))
		data_index = inbuilt_moon.index(closest_X) + 1

		#load sky transmission & extinction files, then reduce to the columns required
		moon_em_all_X = np.genfromtxt(os.path.join(bgpath, 'moon.txt'), delimiter=',')
		
		moon_em_lambda = moon_em_all_X[:,0]
		moon_em_flux = moon_em_all_X[:,data_index]

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
def sky_transmission(lambs, air_mass, debug_plots, output_file):
	'''Function that generates a full throughput curve combining
	sky transmission & sky extinctionp.

	Inputs:
		lambs: array of wavelengths for datacube
		air_mass: Air mass of the observation

	Outputs:
		cube_total_sky_trans: array of total throughput
			for each wavelength value in lambs
	'''
	#convert from zenith angle to airmass
	inbuilt_airmasses = [1.1, 1.3, 1.5, 2.0]
	if air_mass not in inbuilt_airmasses:
		raise ValueError('Error: ' + str(air_mass) + ' is not a valid air_mass. Valid options are: ' + ",".join([str(_) for _ in inbuilt_airmasses]))

	#determine the closest data to the airmass value given and find it's location in the data file
	closest_X = min(inbuilt_airmasses, key=lambda x:abs(x - air_mass))
	data_index = inbuilt_airmasses.index(closest_X) + 1

	#load sky transmission & extinction files, then reduce to the columns required
	sky_trans_all_X = np.genfromtxt(os.path.join(bgpath, 'transmission.txt'), delimiter=',')
	
	sky_tr_lambda = sky_trans_all_X[:,0]
	sky_tr = sky_trans_all_X[:,data_index]

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
	



def sim_sky(cube, back_emission, header, ext_lambs, cube_lamb_mask, DIT, air_mass, moon, site_temp, adr_switch, debug_plots=False, output_file=""):
	''' Simulates sky effects
	Inputs:
		cube: Input datacube (RA, DEC, lambda)
		back_emission: Input background emission outside of the FoV
		header: FITS header
		ext_lambs: extended lambda array [um]
		cube_lamb_mask: mask array to get the lambs of the cube
		DIT: Exposure time [s]
		air_mass: Air mass of the observation
		moon: Fractional moon illumination
		site_temp: Telescope temperature [K]
		adr_switch: Boolean - turn ADR on or off
		debug_plots: Produce debug plots
		output_file: File name for debug plots
	Outputs:
		cube: Cube including sky emission, transmission and ADR
		back_emission: back_emission including sky
	'''
	
	# Get sky transmission
	print "Calculating sky transmission"
	sky_trans = sky_transmission(ext_lambs, air_mass, debug_plots, output_file)
	
	# Get sky emission (lines + continuum)
	print "Calculating sky emission"
	sky_emission = sky_background(ext_lambs, air_mass, DIT, debug_plots, output_file)
	
	# Get moon emission
	print "Calculating Moon emission"
	moon_emission = moon_background(ext_lambs, moon, DIT, debug_plots, output_file)
	back_emission = back_emission + sky_emission + moon_emission
	

	# Add sky emission/transmission to the input cube
	sky_trans_cube = sky_trans[cube_lamb_mask]
	sky_trans_cube.shape = (np.sum(cube_lamb_mask),1,1)
	cube = np.multiply(cube, sky_trans_cube)

	sky_emission_cube = sky_emission[cube_lamb_mask] + moon_emission[cube_lamb_mask]
	sky_emission_cube.shape = (np.sum(cube_lamb_mask),1,1)
	cube = np.add(cube, sky_emission_cube)
	
	# Add atmospheric differential refration
	if adr_switch == "True":
		print "Calculating ADR"
		lambs = ext_lambs[cube_lamb_mask]
		cube = apply_adr(cube, header, lambs, site_temp, air_mass, debug_plots=False, output_file=output_file)
		
		
	return cube, back_emission

