# coding=utf-8
'''
Helper functions to model the emissivity
'''
import os

import numpy as np
import scipy.constants as sp
from scipy.interpolate import interp1d

from src.config import *
from src.modules.misc_utils import path_setup

import matplotlib.pylab as plt

tppath = path_setup('../../' + config_data["data_dir"] + 'throughput/')

def blackbody(waves, T):
	'''Function that gives the Plank blackbody spectrum
	as a function of wavelength and temperature.

	Inputs:
		waves: value or array of wavelengths in microns
		T: Temperature in Kelvin

	Outputs:
		bb_spectrum: value or array for Plank blackbody spectrum
			units - [J/s/m^2/lambda(um)/arcsec^2]
	'''

	#Convert wavelength array from microns to metres
	wave = waves*1.E-6
	
	exp_part = np.exp(sp.h*sp.c/(wave*sp.k*T))

	#Flux in [J/s/m/m2/steradian]
	bb_spectrum = (2.*sp.h*sp.c**2/wave**5)*(exp_part - 1)**(-1)
	#put into units of: J/s/lambda(um)/m^2/arcsec^2
	bb_spectrum /= 1.E6 #to get into J/s/m2/um/steradian
	bb_spectrum /= 4.2545E10 #to get into J/s/m2/um/arcsec2
	

	return bb_spectrum



def load_transmission_curve(wavels, filename, show_plot, plot_file, plot_label):
	'''Load a transmission curve from a file

	Inputs:
		wavels: array of wavelengths for datacube
		filename: file containing the transmission curve
		show_plot: plot curve. True/False
		plot_file: vector with the name of the debug plot
		plot_label: y label of the plot

	Outputs:
		cube_trans: array of throughput for each wavelength value in wavels
	'''

	data = np.genfromtxt(os.path.join(tppath, filename), delimiter=',')

	#Interpolate as a function of wavelength
	trans_interp = interp1d(data[:, 0], data[:, 1],
				kind='linear', bounds_error=False, fill_value=0.)
	#Obtain values for datacube wavelength array
	cube_trans = trans_interp(wavels)
	
	
	if show_plot:
		plt.clf()
		plt.plot(wavels, cube_trans)
		plt.xlabel(r"wavelength [$\mu$m]")
		plt.ylabel(plot_label)
		plt.savefig("_".join(plot_file) + "_tr.pdf")
		np.savetxt("_".join(plot_file) + "_tr.txt", np.c_[wavels, cube_trans])

	return cube_trans


def get_background_emission(wavels, T, emissivity, DIT, show_plot, plot_file, plot_label):
	'''Function that calculates the background emission
	using emissivity and the Plank BB function.
	Emissivity modelled as a graybody with eps*BB(T)
	
	Inputs:
		wavels: array of wavelengths for datacube
		T: temperature [K]
		emissivity: emissivity typically modelled as 1 - reflectivity
		DIT: exposure time [s].
		show_plot: plot curve. True/False
		plot_file: vector with the name of the debug plot
		plot_label: y label of the plot
		
	Outputs:
		bg_spec_ph: array of total background emission
			[units of photons/m^2/um/arcsec^2]
			for each wavelength value in wavels
	'''
	
	cube_bb_spec = blackbody(wavels, T)

	bg_spec = np.multiply(emissivity, cube_bb_spec)

	bg_spec_ph = bg_spec/(sp.h*sp.c/(wavels*1.E-6))*DIT

	if show_plot:
		plt.clf()
		plt.plot(wavels, bg_spec_ph, label="Blackbody T = {:.1f} K".format(T))
		plt.legend()
		plt.xlabel(r"wavelength [$\mu$m]")
		plt.ylabel(plot_label)
		plt.savefig("_".join(plot_file) + "_em.pdf")
		np.savetxt("_".join(plot_file) + "_em.txt", np.c_[wavels, bg_spec_ph])
	
	return bg_spec_ph

