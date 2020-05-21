'''Module to calculate and apply atmospheric differential refraction
to datacubes.
'''

from src.modules.rebin import *
import numpy as np
import matplotlib.pylab as plt
from scipy.interpolate import interp2d

def apply_adr(cube, header, wave, T, am, correct=False, debug_plots=False, output_file=""):
	''' calculates the differential refration in MILLI-ARCSEC
	across the wavelength range provided given the refractive
	index and airmass.

	INPUTS:
		cube: data cube
		header: data cube header
		wave:       wavelength axis, in microns
		T:      temperature (Kelvin)
		am:         airmass, unitless (= sec z)
		correct: if True, correct for ADR in the data cube
		debug_plots: Produce debug plots
		output_file: File name for debug plots

	'''
	
	refr = calc_ref(wave, T)
	optlam = optimalguide(wave[0], wave[-1], T)
	adr = calc_adr(wave, refr, am, optlam)

	if debug_plots:
		plt.clf()
		plt.plot(wave, adr, label="ADR T = {:.1f} K airmass = {:.1f}".format(T, am))
		plt.legend()
		plt.xlabel(r"wavelength [$\mu$m]")
		plt.ylabel(r"ADR [mas]")
		plt.savefig(output_file + "_adr.pdf")

	
	# apply ADR to the input cube
	z, y, x = cube.shape
		
	xmax = header['CDELT1']*(x-1)
	ymax = header['CDELT2']*(y-1)

	xgrid_out = np.linspace(0, xmax, x)
	ygrid_out = np.linspace(0, ymax, y)
	
	for k in np.arange(0, z):
		xgrid_in = xgrid_out
		if correct == False:
			ygrid_in = ygrid_out + adr[k]
		else:
			ygrid_in = ygrid_out - adr[k]
		
		image = interp2d(xgrid_in, ygrid_in, cube[k,:,:], kind='cubic')
		cube[k,:,:] = image(xgrid_out, ygrid_out)

	return cube


def calc_ref(l, T, rh=0.1, P=760.0):
	'''This function calculates the refractive index of air vs. wavelength.
	ppwv derived from relative humidity RH, saturation
	pressure of water ps, and temperature T. Default values for RH and P
	from Kendrew et al. 2008.
	EQUATIONS:

	RH = ppwv / ps
	ps = 6.11 * exp([17.27 * T] / [T + 237.3]) [millibar] for T in deg C.
	Reference: VDF-TRE-IOA-00009-0003 by Dafydd Wyn Evans, 2004, http://casu.ast.cam.ac.uk/documents/wfcam/astrometry/refractir.pdf

	INPUTS:
	l:      wavelength axis (microns)
	T:      temperature (Kelvin)
	rh:     relative humidity (fraction)
	P:      atmospheric pressure (millibar)

	OUTPUT:
	refractive index of air on the provided wavelength scale'''

	Tc = T-273.15
	ppwv = rh * 6.11 * np.exp((17.27 * Tc) / (Tc + 237.3))

	# Ps = 1013.25
	Ps = 1000.0#When using mbar units --> Reference: VDF-TRE-IOA-00009-0003 by Dafydd Wyn Evans, 2004
	Ts = 288.15

	nout = 1 + ( 64.328 + (29498.1)/(146.0 - l**(-2)) + (255.4)/(41 - l**(-2)) )\
		* ((P*Ts)/(Ps*T)) * 1e-6 - 43.49 * (1.0 - (7.956e-3)/l**2)*(ppwv/Ps) * 1e-6

	return nout

def calc_adr(wave, n, am, optlam):

	''' calculates the differential refration in MILLI-ARCSEC
	across the wavelength range provided given the refractive
	index and airmass.

	INPUTS:
	wave:       wavelength axis, in microns
	n:          refractive index, unitless (must have same length as wave)
	am:         airmass, unitless (= sec z)
	optlam:     guiding wavelength, in microns
	'''

	wavearg = np.where(wave > optlam)[0][0]
	nc = n[wavearg]
	#print 'Guiding wavelegth: %.3f' % optlam

	z = np.arccos(1./am)

	diffref = 206265. * 1.e3 * ((n**2-1.0)/(2.0*n**2) - (nc**2-1.0)/(2.*nc**2)) * np.tan(z)

	return diffref



def optimalguide(wave0, wave1, temp) :
	'''
	Calculates the optimal guiding wavelength for to given end wavelengths
	'''

	# Only need to work with refractive index, as that's the only wavelength dependant term

	nLow = calc_ref(wave0, temp)
	nHigh = calc_ref(wave1, temp)
	nMid = nLow + (nHigh-nLow)/2.0

	# iterate to find the optimal wavelength

	diff=0.0
	dw = wave1-wave0
	waveg = wave0 + (wave1-wave0)/2.0

	while (True) : # 1% tolerance on finding where nMid is??
		waveg=min(max(wave0,waveg + min(1.0,dw)*diff*100000.0),wave1)
		nwave = calc_ref(waveg, temp)
		diff=nwave-nMid
		if(np.abs(diff) < (np.abs(nHigh-nLow))*0.001) : break

	#print 'Guiding wavelength = ', waveg
	return waveg
