'''
Code for the bulk framework of the HARMONI
simulator. This code should show the main functions/processes
to move from an input datacube (lambda, y, x) to output cubes:
'''

import collections
import datetime
import multiprocessing as mp
import os.path

import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import convolve2d
from astropy.io import fits

import matplotlib.pylab as plt

from config import *
from init_cube import init_cube
from sim_sky import sim_sky
from sim_telescope import sim_telescope
from sim_instrument import sim_instrument
from sim_detector import sim_detector, apply_crosstalk
from modules.adr import apply_adr

from modules.rebin import *

def main(datacube, outdir, DIT, NDIT, grating, spax, seeing, air_mass, version, res_jitter=3., moon=0.,
	 site_temp=280.5, adr_switch='True', seednum=100, nprocs=mp.cpu_count()-1, keep_debug_plots=False, aoMode="LTAO"):
	'''
	Inputs:
		datacube: Input high resolution datacube (RA, DEC, lambda)
		outdir: Output file directory
		DIT: Exposure time [s]
		NDIT: No. of exposures
		grating: Spectral grating
		spax: spatial pixel (spaxel) scale [mas]
		seeing: Atmospheric seeing FWHM [arcsec]
		air_mass: Air mass of the observation
		res_jitter: Residual telescope jitter [mas]
		site_temp: Telescope temperature [K]
		moon: Fractional Moon illumination
		adr_switch: Boolean - turn ADR on or off.
		seednum: ramdom seed number
		nprocs: Number of processes
		keep_debug_plots: keep debug plots
		aoMode: Adaptive optics mode: "LTAO" or "noAO" for seeing limited

	Outputs:

	'''
	debug_plots = True
	#keep_debug_plots = True
	
	Conf = collections.namedtuple('Conf', 'name, header, value')
	simulation_conf = [
			Conf('HSIM Version', 'HSM_VER', version),
			Conf('Filename', 'HSM_FILE', datacube),
			Conf('Output dir', 'HSM_OUTD', outdir),
			Conf('DIT', 'HSM_DIT', DIT),
			Conf('NINT', 'HSM_NINT', NDIT),
			Conf('Spaxel', 'HSM_SPAX', spax),
			Conf('Grating', 'HSM_GRAT', grating),
			Conf('Seeing', 'HSM_SEEI', seeing),
			Conf('Air Mass', 'HSM_AIRM', air_mass),
			Conf('Residual telescope jitter', 'HSM_JITT', res_jitter),
			Conf('Temperature', 'HSM_TEMP', site_temp),
			Conf('Moon', 'HSM_MOON', moon),
			Conf('ADR', 'HSM_ADR', adr_switch),
			Conf('Seed', 'HSM_SEED', seednum),
			Conf('AO', 'HSM_AO', aoMode),
			Conf('No. of processes', 'HSM_NPRC', nprocs),
			]
	
	for _ in simulation_conf:
		print _.name + " = " + str(_.value)
	print 
	
	if aoMode.upper() == "LTAO":
		use_LTAO = True
	elif aoMode.upper() == "NOAO":
		use_LTAO = False
	else:
		raise ValueError('Error: ' + aoMode + ' is not a valid AO mode. Valid options are: LTAO or noAO')

	
	base_name = os.path.splitext(os.path.basename(datacube))[0]
	base_filename = os.path.join(outdir, base_name)

	np.random.seed(seednum)

	# Read input FITS cube and resample depending on grating and spaxel scale
	# output is in ph/s/m2/um/arcsec2 units
	cube, head, lambs = init_cube(datacube, grating, spax)
	
	# Calculate extended lambda range to convolve with the LSF
	LSF_width = int(config_data["spectral_sampling"]["internal"]/2.35482*config_data['LSF_kernel_size'])
	lambs_extended_index = (np.arange(LSF_width)+1)*(lambs[1] - lambs[0])
	
	lambs_extended = np.concatenate(((lambs[0] - lambs_extended_index)[::-1], lambs, lambs[-1] + lambs_extended_index))
	cube_lamb_mask = np.concatenate(([False]*len(lambs_extended_index), [True]*len(lambs), [False]*len(lambs_extended_index)))
	
	# calculate the cube in photons for a single exposure
	# in ph/m2/um/arcsec2
	cube_exp = cube*DIT
	back_emission = np.zeros(len(lambs_extended))
	
	# 1 - Atmosphere: 
	#	- Sky emission (lines + continuum)
	#	- Sky transmission
	#	- Moon
	#	- Atmospheric differential refration
	cube_exp, back_emission = sim_sky(cube_exp, back_emission, head, lambs_extended, cube_lamb_mask, DIT, air_mass, moon, site_temp, adr_switch, debug_plots=debug_plots, output_file=base_filename)
	
	# 2 - Telescope:
	#	- PSF + Jitter
	#	- Telescope background
	#	- Telescope transmission
	cube_exp, back_emission, psf = sim_telescope(cube_exp, back_emission, lambs_extended, cube_lamb_mask, DIT, res_jitter, air_mass, seeing, spax, site_temp, nprocs, debug_plots=debug_plots, output_file=base_filename, use_LTAO=use_LTAO)

	# 3 - Instrument:
	#	- LSF
	#	- Instrument background
	#	- Instrument transmission
	cube_exp, back_emission = sim_instrument(cube_exp, back_emission, lambs_extended, cube_lamb_mask, DIT, grating, site_temp, debug_plots=debug_plots, output_file=base_filename, LTAO_dichroic=use_LTAO)
		
	# 4 - Rebin cube to output spatial and spectral pixel size
	print "Rebin data"
	print "- Output spaxel scale: ", spax 
	z, y, x = cube_exp.shape
	
	# rebin spatial axes
	spax_scale = config_data['spaxel_scale'][spax]
	scale_x = spax_scale.xscale/spax_scale.psfscale
	scale_y = spax_scale.yscale/spax_scale.psfscale
	out_size_x = int(x/scale_x)
	out_size_y = int(y/scale_y)
	output_cube = np.zeros((z, out_size_y, out_size_x))
	
	for k in np.arange(0, z):
		output_cube[k,:,:] = frebin2d(cube_exp[k,:,:], (out_size_x, out_size_y))
	
	# and update header
	head['CDELT1'] = spax_scale.xscale
	head['CDELT2'] = spax_scale.yscale
		
	# rebin spectral axis
	scale_z = config_data["spectral_sampling"]["output"]/config_data["spectral_sampling"]["internal"]
	
	lambs = lambs_extended[cube_lamb_mask]
	new_lamb_per_pix = (lambs[1] - lambs[0])/scale_z
	
	output_lambs = new_lamb_per_pix*np.arange(int(len(lambs)*scale_z)) + lambs[0]
	output_cube_spec = np.zeros((len(output_lambs), out_size_y, out_size_x))
	
	print "- Output spectral sampling: {:.2f} A".format(new_lamb_per_pix*10000.)
	
	for i in np.arange(0, out_size_x):
		for j in np.arange(0, out_size_y):
			output_cube_spec[:,j,i] = rebin1d(output_lambs, lambs, output_cube[:,j,i])
	
	output_back_emission = rebin1d(output_lambs, lambs_extended, back_emission)
	
	# and update header
	head['CRVAL3'] = output_lambs[0]
	head['CDELT3'] = new_lamb_per_pix
	head['NAXIS3'] = len(output_lambs)
	
	
	# correct for ADR
	if adr_switch == "True":
		print "- Correcting ADR"
		output_cube_spec = apply_adr(output_cube_spec, head, output_lambs, site_temp, air_mass, correct=True)

	
	# 5 - Convert to photons per pixel
	spaxel_area = spax_scale.xscale/1000.*spax_scale.yscale/1000. # [arcsec2]
	channel_width = new_lamb_per_pix
	
	output_cube_spec = output_cube_spec*spaxel_area*channel_width*config_data["telescope"]["area"]
	output_back_emission = output_back_emission*spaxel_area*channel_width*config_data["telescope"]["area"]
	head['FUNITS'] = "photons"
	
	# 6 - Detector
	#	- QE
	#	- Dark
	#	- Read noise
	output_cube_spec, output_back_emission, read_noise, dark_current = sim_detector(output_cube_spec, output_back_emission, output_lambs, grating, DIT, debug_plots=debug_plots, output_file=base_filename)
	head['FUNITS'] = "electrons"
	
	
	# Generate noiseless outputs
	# - object + background
	output_cube_spec = output_cube_spec.astype(np.float32)
	
	# - background - as a cube the 1D background spectrum
	output_back_emission = output_back_emission.astype(np.float32)
	output_back_emission.shape = (len(output_back_emission),1,1)
	output_back_emission_cube = np.zeros_like(output_cube_spec) + output_back_emission
	
	# - object - back
	output_cube_spec_wo_back = np.subtract(output_cube_spec, output_back_emission_cube)

	# Generate observed outputs
	zero_cube = np.zeros_like(output_cube_spec)
	
	# A. Object exposure
	# - object exposure with crosstalk
	sim_object_plus_back = np.random.poisson(abs(output_cube_spec*NDIT)).astype(np.float32)
	sim_object_plus_back = apply_crosstalk(sim_object_plus_back, config_data["crosstalk"])

	# - read noise and dark current for object exposure
	dark_cube = np.zeros_like(output_cube_spec) + dark_current*NDIT
	sim_read_noise1 = np.random.normal(zero_cube, np.sqrt(NDIT)*read_noise).astype(np.float32)
	sim_dark_current1 = np.random.poisson(dark_cube).astype(np.float32)
	
	# - combine object, read noise and dark current
	sim_total = sim_object_plus_back + sim_read_noise1 + sim_dark_current1
	
	
	# B. Background exposure
	#- background. No crosstalk since the background is uniform
	sim_back = np.random.poisson(abs(output_back_emission_cube*NDIT)).astype(np.float32)
	
	# - read noise and dark current for background exposure
	sim_read_noise2 = np.random.normal(zero_cube, np.sqrt(NDIT)*read_noise).astype(np.float32)
	sim_dark_current2 = np.random.poisson(dark_cube).astype(np.float32)
	
	# - combine object, read noise and dark current
	sim_total_only_back = sim_back + sim_read_noise2 + sim_dark_current2
	
	
	# C. Calculate reduced cube: object - background exposure
	sim_reduced = sim_total - sim_total_only_back
	
	
	
	print "Pipeline interpolation effects"
	# Convolve the reduced cube with a 1pix FWHM Gaussian to account for the 
	# interpolation during the data reduction
	sigma = 1./2.35482 # pix
	kernel_size = 5
	Gauss2D = lambda x, y: 1./(2.*np.pi*sigma**2)*np.exp(-(x**2 + y**2)/(2.*sigma**2))
	xgrid = np.linspace(1, kernel_size, kernel_size) - kernel_size*0.5 - 0.5
	ygrid = np.linspace(1, kernel_size, kernel_size) - kernel_size*0.5 - 0.5
	xx, yy = np.meshgrid(xgrid, ygrid)
	kernel = Gauss2D(xx, yy)
	z, y, x = sim_reduced.shape
	for pz in range(z):
		sim_reduced[pz,:,:] = convolve2d(sim_reduced[pz,:,:], kernel, mode="same", boundary="wrap")

	
	# Calculate noise cube
	noise_cube_object = abs(output_cube_spec*NDIT) # object+back noise variance
	noise_cube_back = abs(output_back_emission_cube*NDIT) # back noise variance
	noise_cube_read_noise = np.sqrt(NDIT)*read_noise # read noise sigma
	noise_cube_dark = dark_cube # dark noise variance
	
	noise_cube_total = np.sqrt(noise_cube_object + noise_cube_back + 2.*noise_cube_dark + 2.*noise_cube_read_noise**2)
	

	#
	print "Saving output"
	if debug_plots:
		plt.clf()

		colors = [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf']

		fig = plt.figure()
		ax = fig.add_subplot(111)
		
		total_telescope_sky_em = np.zeros_like(lambs_extended)
		total_instrument_em = np.zeros_like(lambs_extended)
		
		w, e = np.loadtxt(base_filename + "_sky_em.txt", unpack=True)
		plt.plot(w, e, label="sky", color=colors[0])
		total_telescope_sky_em += e
		w, e = np.loadtxt(base_filename + "_tel_em.txt", unpack=True)
		plt.plot(w, e, label="telescope", color=colors[1])
		total_telescope_sky_em += e

		if moon > 0.:
			w, e = np.loadtxt(base_filename + "_moon_em.txt", unpack=True)
			plt.plot(w, e, label="Moon", color=colors[6])
			total_telescope_sky_em += e
		
		if use_LTAO:
			w, e = np.loadtxt(base_filename + "_ins_LTAOd_em.txt", unpack=True)
			plt.plot(w, e, label="LTAO dichroic", color=colors[2])
			total_instrument_em += e
		
		w, e = np.loadtxt(base_filename + "_ins_FPRS_em.txt", unpack=True)
		plt.plot(w, e, label="FPRS", color=colors[3])
		total_instrument_em += e
		w, e = np.loadtxt(base_filename + "_ins_em.txt", unpack=True)
		plt.plot(w, e, label="instrument after FPRS", color=colors[4])
		total_instrument_em += e
		
		plt.legend()
		plt.xlabel(r"wavelength [$\mu$m]")
		plt.ylabel(r"back emission [photons/m$^2$/$\mu$m/arcsec$^2$]")
		plt.yscale("log")
		plt.text(0.1, 0.2, "HARMONI/(Telescope+Sky) = {:.2f}".format(np.nanmedian(total_instrument_em/total_telescope_sky_em)), transform = ax.transAxes)
		plt.savefig(base_filename + "_total_em.pdf")

		plt.clf()
		w, e = np.loadtxt(base_filename + "_sky_tr.txt", unpack=True)
		plt.plot(w, e, label="sky", color=colors[0])
		w, e = np.loadtxt(base_filename + "_tel_tr.txt", unpack=True)
		plt.plot(w, e, label="telescope", color=colors[1])
		if use_LTAO:
			w, e = np.loadtxt(base_filename + "_ins_LTAOd_tr.txt", unpack=True)
			plt.plot(w, e, label="LTAO dichroic", color=colors[2])
		w, e = np.loadtxt(base_filename + "_ins_FPRS_tr.txt", unpack=True)
		plt.plot(w, e, label="FPRS", color=colors[3])
		w, e = np.loadtxt(base_filename + "_ins_tr.txt", unpack=True)
		plt.plot(w, e, label="instrument after FPRS", color=colors[4])
		w, e = np.loadtxt(base_filename + "_det_qe.txt", unpack=True)
		plt.plot(w, e, label="detector", color=colors[5])
		
		plt.legend()
		plt.xlabel(r"wavelength [$\mu$m]")
		plt.ylabel(r"transmission")
		plt.savefig(base_filename + "_total_tr.pdf")

		if not keep_debug_plots:
			list_files = ["sky_tr", "sky_em", "moon_em", "tel_tr", "tel_em", "ins_tr", "ins_em", "det_qe", "ins_FPRS_tr", "ins_FPRS_em"]
			if use_LTAO:
				list_files.append("ins_LTAOd_tr")
				list_files.append("ins_LTAOd_em")
			for _ in list_files:
				os.remove(base_filename + "_" + _ + ".txt")
				os.remove(base_filename + "_" + _ + ".pdf")
		
		
	# Noiseless
	outFile_noiseless_object = base_filename + "_noiseless_obj.fits"
	outFile_noiseless_background = base_filename + "_noiseless_back.fits"
	outFile_noiseless_object_plus_back = base_filename + "_noiseless_obj_plus_back.fits"
	# With noise.
	# - Observed Obj + Back + Noise1
	# - Background  Back + Noise2
	# - Reduced  Obj + Back + Noise1 - (Back + Noise2)
	outFile_observed = base_filename + "_observed_obj_plus_back.fits"
	outFile_observed_back = base_filename + "_observed_back.fits"
	outFile_reduced = base_filename + "_reduced.fits"
	# SNR
	# - Obj/Noise
	outFile_SNR = base_filename + "_reduced_SNR.fits"
	
	# Update header
	for _ in simulation_conf:
		head[_.header] = str(_.value)
	
	head['HSM_TIME'] = str(datetime.datetime.utcnow())
	
	save_fits_cube(outFile_noiseless_object_plus_back, output_cube_spec*NDIT, "Noiseless O+B", head)
	save_fits_cube(outFile_noiseless_background, output_back_emission_cube*NDIT, "Noiseless B", head)
	save_fits_cube(outFile_noiseless_object, output_cube_spec_wo_back*NDIT, "Noiseless O", head)
	
	save_fits_cube(outFile_observed, sim_total, "Observed O+B+Noise", head)
	save_fits_cube(outFile_reduced, sim_reduced, "Reduced (O+B+Noise) - (B+Noise)", head)
	
	save_fits_cube(outFile_SNR, output_cube_spec_wo_back*NDIT/noise_cube_total, "SNR (O-B)/Noise", head)
	
	# Save PSF
	head_PSF = head
	head_PSF['CDELT1'] = spax_scale.psfscale
	head_PSF['CDELT2'] = spax_scale.psfscale
	head_PSF['LAMBDA'] = lambs_extended[cube_lamb_mask][0]
	try:
		del head_PSF['CDELT3']
		del head_PSF['CRVAL3']
		del head_PSF['CRPIX3']
		del head_PSF['CTYPE3']
		del head_PSF['CUNIT3']
		del head_PSF['FUNITS']
		del head_PSF['SPECRES']
	except:
		pass
	
	save_fits_cube(base_filename + "_PSF.fits", psf, "PSF", head_PSF)


	print 'Simulation Complete'
	return
	

def save_fits_cube(filename, data, typ, header):
	header['HSM_TYPE'] = typ
	fits.writeto(filename, data, header=header, overwrite=True)
	
	
