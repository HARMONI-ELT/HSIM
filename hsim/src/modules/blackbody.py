import numpy as np
import scipy.constants as sp

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