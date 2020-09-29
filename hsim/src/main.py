'''
Code for the bulk framework of the HARMONI
simulator. This code should show the main functions/processes
to move from an input datacube (lambda, y, x) to output cubes:
'''
import collections
import glob
import re
import datetime
import multiprocessing as mp
import os.path
import logging
import warnings

import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import convolve2d
from astropy.io import fits
import astropy.constants as const
from photutils import aperture_photometry, CircularAperture

import matplotlib.pylab as plt

from src.config import *
from src.init_cube import init_cube
from src.sim_sky import sim_sky
from src.sim_telescope import sim_telescope
from src.sim_instrument import sim_instrument
from src.sim_detector import sim_detector, apply_crosstalk, mask_saturated_pixels, apply_crosstalk_1d
from src.modules.adr import apply_adr
from src.modules.rebin import *

# Detector systematics
from src.sim_detector import make_rn_dist, make_dets, add_detectors
from src.modules.misc_utils import trim_cube


def main(datacube, outdir, DIT, NDIT, grating, spax, seeing, air_mass, version,\
		 res_jitter=3., moon=0., site_temp=280.5, adr_switch='True', \
		 det_switch='False', det_save_path="None", seednum=100, nprocs=mp.cpu_count()-1, \
		 debug=True, aoMode="LTAO"):
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
		res_jitter: Residual telescope jitter. 1 or 2 x separated numbers [mas]
		site_temp: Telescope temperature [K]
		moon: Fractional Moon illumination
		adr_switch: Boolean - turn ADR on or off.
		det_switch: Boolean - use detector systematics (off by default)
		det_save_path: Directory to save interim detector files
		seednum: ramdom seed number
		nprocs: Number of processes
		debug: keep debug plots and all noise outputs
		aoMode: Adaptive optics mode: "LTAO", "SCAO", "Airy" or "noAO" for seeing limited

	Outputs:

	'''
	warnings.filterwarnings("ignore", module="astropy")
	warnings.filterwarnings("ignore", module="matplotlib")
	debug_plots = True

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
			Conf('Detectors', 'HSM_DET', det_switch),
			Conf('Seed', 'HSM_SEED', seednum),
			Conf('AO', 'HSM_AO', aoMode.upper()),
			Conf('No. of processes', 'HSM_NPRC', nprocs),
			]

	base_name = os.path.splitext(os.path.basename(datacube))[0]
	base_filename = os.path.join(outdir, base_name)

	logfile = base_filename + ".log"
	open(logfile, 'w').close()

	# log to a file
	logging.basicConfig(filename=logfile, level=logging.INFO, format='%(asctime)s  %(levelname)s  %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
	# and also to the terminal
	logger = logging.getLogger()
	std = logging.StreamHandler()
	std.setFormatter(HSIMFormatter())
	std.setLevel(logging.INFO)
	logger.addHandler(std)

	hsimlog = HSIMLoggingHandler()
	logger.addHandler(hsimlog)


	logging.info("Simulation input parameters:")
	for _ in simulation_conf:
		logging.info(_.name + " = " + str(_.value))

	if aoMode.upper() in ["LTAO", "SCAO", "NOAO", "AIRY"]:
		aoMode = aoMode.upper()
	elif not os.path.isfile(aoMode):
		logging.error(aoMode + ' is not a valid AO mode. Valid options are: LTAO, SCAO, noAO, Airy')
		return

	if air_mass not in config_data["PSD_cube"]["air_masses"]:
		logging.error(str(air_mass) + ' is not a valid air mass. Valid options are: ' + ", ".join(map(str, sorted(config_data["PSD_cube"]["air_masses"]))))
		return

	if seeing not in  config_data["PSD_cube"]["seeings"]:
		logging.error(str(seeing) + ' is not a valid seeing. Valid options are: ' + ", ".join(map(str, sorted(config_data["PSD_cube"]["seeings"]))))
		return


	try:
		if "x" in res_jitter:
			jitter = np.array([*map(float, res_jitter.split("x"))])
			if len(jitter) != 2:
				logging.error(str(res_jitter) + " is not a valid jitter value.")
				return
		else:
			jitter = np.repeat(float(res_jitter), 2)
		
	except:
		logging.error(str(res_jitter) + " is not a valid jitter value.")
		return

	np.random.seed(seednum)

	# Read input FITS cube and resample depending on grating and spaxel scale
	# output is in ph/s/m2/um/arcsec2 units
	cube, head, lambs, input_spec_res = init_cube(datacube, grating, spax)

	# Calculate extended lambda range to convolve with the LSF
	LSF_width = int(config_data["spectral_sampling"]["internal"]/2.35482*config_data['LSF_kernel_size'])
	lambs_extended_index = (np.arange(LSF_width)+1)*(lambs[1] - lambs[0])

	lambs_extended = np.concatenate(((lambs[0] - lambs_extended_index)[::-1], lambs, lambs[-1] + lambs_extended_index))
	cube_lamb_mask = np.concatenate(([False]*len(lambs_extended_index), [True]*len(lambs), [False]*len(lambs_extended_index)))

	# calculate the cube in photons for a single exposure
	# in ph/m2/um/arcsec2
	cube_exp = cube*DIT
	back_emission = np.zeros(len(lambs_extended)) # Background exposure
	transmission = np.ones(len(lambs_extended)) # Telluric correction

	# 1 - Atmosphere: 
	#	- Sky emission (lines + continuum)
	#	- Sky transmission
	#	- Moon
	#	- Atmospheric differential refration
	cube_exp, back_emission, transmission = sim_sky(cube_exp, back_emission, transmission, head, lambs_extended, cube_lamb_mask, DIT, air_mass, moon, site_temp, adr_switch, input_spec_res, debug_plots=debug_plots, output_file=base_filename)
	
	# 2 - Telescope:
	#	- PSF + Jitter
	#	- Telescope background
	#	- Telescope transmission
	cube_exp, back_emission, transmission, psf_internal, psf_lambda = sim_telescope(cube_exp, back_emission, transmission, lambs_extended, cube_lamb_mask, DIT, jitter, air_mass, seeing, spax, site_temp, aoMode, nprocs, debug_plots=debug_plots, output_file=base_filename)

	# 3 - Instrument:
	#	- LSF
	#	- Instrument background
	#	- Instrument transmission
	cube_exp, back_emission, transmission, LSF_width_A = sim_instrument(cube_exp, back_emission, transmission, lambs_extended, cube_lamb_mask, DIT, grating, site_temp, input_spec_res, aoMode, debug_plots=debug_plots, output_file=base_filename)
	
	# 4 - Rebin cube to output spatial and spectral pixel size
	logging.info("Rebin data")
	logging.info("- Output spaxel scale: " + str(spax))
	z, y, x = cube_exp.shape
	
	# rebin spatial axes
	spax_scale = config_data['spaxel_scale'][spax]
	scale_x = spax_scale.xscale/spax_scale.psfscale
	scale_y = spax_scale.yscale/spax_scale.psfscale
	out_size_x = int(x/scale_x)
	out_size_y = int(y/scale_y)
	output_cube = np.zeros((z, out_size_y, out_size_x))
	
	for k in np.arange(0, z):
		output_cube[k, :, :] = frebin2d(cube_exp[k, :, :], (out_size_x, out_size_y))
	
	# and update header
	head['CDELT1'] = spax_scale.xscale
	head['CDELT2'] = spax_scale.yscale
		
	# rebin spectral axis
	scale_z = config_data["spectral_sampling"]["output"]/config_data["spectral_sampling"]["internal"]
	
	lambs = lambs_extended[cube_lamb_mask]
	new_lamb_per_pix = (lambs[1] - lambs[0])/scale_z # micron
	
	output_lambs = new_lamb_per_pix*np.arange(int(len(lambs)*scale_z)) + lambs[0]
	
	logging.info("- Output spectral sampling: {:.2f} AA".format(new_lamb_per_pix*10000.))
	
	output_cube_spec = rebin_cube_1d(output_lambs, lambs, output_cube)
	output_back_emission = rebin1d(output_lambs, lambs_extended, back_emission)
	output_transmission = rebin1d(output_lambs, lambs_extended, transmission)
	
	# and update header
	head['CRVAL3'] = output_lambs[0]
	head['CDELT3'] = new_lamb_per_pix
	head['NAXIS3'] = len(output_lambs)
	
	
	# correct for ADR
	if adr_switch == "True":
		logging.info("- Correcting ADR")
		output_cube_spec = apply_adr(output_cube_spec, head, output_lambs, site_temp, air_mass, correct=True)

	
	# 5 - Convert to photons per pixel
	spaxel_area = spax_scale.xscale/1000.*spax_scale.yscale/1000. # arcsec2
	channel_width = new_lamb_per_pix # micron
	
	output_cube_spec = output_cube_spec*spaxel_area*channel_width*config_data["telescope"]["area"]
	output_back_emission = output_back_emission*spaxel_area*channel_width*config_data["telescope"]["area"]
	head['BUNIT'] = "photon"
	
	# 6 - Detector
	#	- QE
	#	- Dark
	#	- Read noise
	#	- Thermal background

	# Cut cubes to correct size if using detector systematics and generate detectors 
	if det_switch == "True" and grating != "V+R":
		logging.info("Trimming datacubes to correct size")
		output_cube_spec = trim_cube(output_cube_spec, verbose=True)
	elif det_switch == "True" and grating == "V+R":
		logging.warning("IR detector systematics selected for visibile grating. Ignoring detector systematics.")
		det_switch = "False"
	
	output_cube_spec, output_back_emission, output_transmission, read_noise, \
					dark_current, thermal_background = sim_detector(\
					output_cube_spec, output_back_emission, \
					output_transmission, output_lambs, grating, DIT, \
					debug_plots=debug_plots, output_file=base_filename)
	head['BUNIT'] = "electron"
	
	
	# Generate noiseless outputs
	# - object + background
	output_cube_spec = output_cube_spec.astype(np.float32)
	# - mask saturated pixels
	output_cube_spec, saturated_obj_back = mask_saturated_pixels(output_cube_spec, grating)
	
	# - background - as a cube the 1D background spectrum
	output_back_emission = output_back_emission.astype(np.float32)
	output_back_emission.shape = (len(output_back_emission), 1, 1)
	output_back_emission_cube = np.zeros_like(output_cube_spec) + output_back_emission
	if det_switch == "True":
		output_back_emission_cube = trim_cube(output_back_emission_cube)
	
	# - mask saturated pixels
	output_back_emission_cube, saturated_back = mask_saturated_pixels(output_back_emission_cube, grating)
	
	# - object - back
	output_cube_spec_wo_back = np.subtract(output_cube_spec, output_back_emission_cube)

	# Generate observed outputs
	zero_cube = np.zeros_like(output_cube_spec)
	
	# A. Object exposure
	# - object exposure with crosstalk
	sim_object_plus_back = np.random.poisson(abs(output_cube_spec*NDIT)).astype(np.float32)
	# Apply crosstalk only to NIR detectors
	if grating != "V+R":
		logging.info("Applying detector crosstalk")
		sim_object_plus_back = apply_crosstalk(sim_object_plus_back, config_data["crosstalk"])
	
	if np.sum(saturated_obj_back) > 0:
		logging.warning(str(np.sum(saturated_obj_back)) + " pixels are saturated in the obj + back frames")
		sim_object_plus_back[saturated_obj_back] = np.nan

	dark_cube = np.zeros_like(output_cube_spec) + dark_current*NDIT
	thermal_cube = np.zeros_like(output_cube_spec) + thermal_background*NDIT
	
	# - read noise and dark current for object exposure
	if det_switch == "False":
		sim_read_noise1 = np.random.normal(zero_cube, np.sqrt(NDIT)*read_noise).astype(np.float32)
	else:
		logging.info("Starting advanced detector systematics")
		rn_dist = make_rn_dist(det_save_path)
		logging.info("- adding systematic effects into observation")
		sim_det_systematics1 = make_dets(rn_dist, DIT)[0]*np.sqrt(NDIT)
	sim_dark_current1 = np.random.poisson(dark_cube).astype(np.float32)
	sim_thermal1 = np.random.poisson(thermal_cube).astype(np.float32)
	
	# - combine object, read noise and dark current
	if det_switch == "False":
		sim_total = sim_object_plus_back + sim_read_noise1 + sim_dark_current1 + sim_thermal1
	else:
		logging.info("- adding detectors into datacube")
		sim_object_plus_dets = add_detectors(sim_object_plus_back, sim_det_systematics1)
		sim_total = sim_object_plus_dets + sim_dark_current1 + sim_thermal1
	
	# B. Background exposure
	#- background with crosstalk
	sim_back = np.random.poisson(abs(output_back_emission_cube*NDIT)).astype(np.float32)
	# Apply crosstalk only to NIR detectors
	if grating != "V+R":
		sim_back = apply_crosstalk(sim_back, config_data["crosstalk"])

	if np.sum(saturated_back) > 0:
		logging.warning(str(np.sum(saturated_back)) + " pixels are saturated in the back frames")
		sim_back[saturated_back] = np.nan

	# - read noise and dark current for background exposure
	if det_switch == "False":
		sim_read_noise2 = np.random.normal(zero_cube, np.sqrt(NDIT)*read_noise).astype(np.float32)
	else:
		logging.info("- creating background exposure")
		sim_det_systematics2 = make_dets(rn_dist, DIT)[0]*np.sqrt(NDIT)
	sim_dark_current2 = np.random.poisson(dark_cube).astype(np.float32)
	sim_thermal2 = np.random.poisson(thermal_cube).astype(np.float32)

	# - combine object, read noise and dark current
	if det_switch == "False":
		sim_total_only_back = sim_back + sim_read_noise2 + sim_dark_current2 + sim_thermal2
	else:
		sim_total_only_back = add_detectors(sim_back, sim_det_systematics2) + \
							  sim_dark_current2 + sim_thermal2

		# - create cube of detector noise
		sim_only_dets = np.zeros_like(sim_total)
		sim_only_dets = add_detectors(sim_only_dets, sim_det_systematics1)
		sim_back_dets = np.zeros_like(sim_total)
		sim_back_dets = add_detectors(sim_back_dets, sim_det_systematics2)
		logging.info("- advanced detector systematics complete")
	
	
	# C. Calculate reduced cube: object - background exposure
	sim_reduced = sim_total - sim_total_only_back
	
	
	logging.info("Pipeline interpolation effects")
	# Convolve the reduced cube with a 1pix FWHM Gaussian to account for the 
	# interpolation during the data reduction
	sigma = 1./2.35482 # pix
	kernel_size = 5
	Gauss2D = lambda x, y: np.exp(-(x**2 + y**2)/(2.*sigma**2))
	xgrid = np.linspace(1, kernel_size, kernel_size) - kernel_size*0.5 - 0.5
	ygrid = np.linspace(1, kernel_size, kernel_size) - kernel_size*0.5 - 0.5
	xx, yy = np.meshgrid(xgrid, ygrid)
	kernel = Gauss2D(xx, yy)
	kernel = kernel/np.sum(kernel)
	z, y, x = sim_reduced.shape
	for pz in range(z):
		sim_reduced[pz, :, :] = convolve2d(sim_reduced[pz, :, :], kernel, mode="same", boundary="fill", fillvalue=0.)

	# Calculate noise cube
	noise_cube_object = abs(output_cube_spec*NDIT) # object+back noise variance
	noise_cube_back = abs(output_back_emission_cube*NDIT) # back noise variance
	noise_cube_read_noise = zero_cube + np.sqrt(NDIT)*read_noise # read noise sigma
	noise_cube_dark = dark_cube # dark noise variance
	noise_cube_thermal = thermal_cube # thermal noise variance
	if det_switch == "True":
		if DIT > 120:
			noise_cube_read_noise = zero_cube + np.sqrt(NDIT)*config_data['systematics']['rd']
		else:
			noise_cube_read_noise = zero_cube + np.sqrt(NDIT)*config_data['systematics']['rd_lowexp']
		noise_cube_pedestal = zero_cube + np.sqrt(NDIT)*config_data['systematics']['pedestal']
		noise_cube_c_pink = zero_cube + np.sqrt(NDIT)*config_data['systematics']['c_pink']
		noise_cube_u_pink = zero_cube + np.sqrt(NDIT)*config_data['systematics']['u_pink']
		noise_cube_acn = zero_cube + np.sqrt(NDIT)*config_data['systematics']['acn']
		noise_cube_pca0 = zero_cube + np.sqrt(NDIT)*config_data['systematics']['pca0_amp']
		noise_cube_dets = zero_cube + np.sqrt(NDIT)*sim_only_dets + np.sqrt(NDIT)*sim_back_dets
	
	if grating != "V+R":
		n_observations = 2
	else:
		n_observations = 1 # no dedicated sky observation
	
	if det_switch == "False":    
		noise_cube_total = np.sqrt(noise_cube_object + n_observations*noise_cube_back + n_observations*noise_cube_dark + n_observations*noise_cube_thermal + n_observations*noise_cube_read_noise**2)
	else:
		noise_cube_total = np.sqrt(noise_cube_object + n_observations*noise_cube_back + n_observations*noise_cube_dark + n_observations*noise_cube_thermal + n_observations*noise_cube_read_noise**2 + \
                                   n_observations*noise_cube_pedestal**2 + n_observations*noise_cube_c_pink**2 + n_observations*noise_cube_u_pink**2 + n_observations*noise_cube_acn**2 + \
                                   n_observations*noise_cube_pca0**2)
		noise_cube_total_with_dets = np.sqrt(noise_cube_object + n_observations*noise_cube_back + n_observations*noise_cube_dark + n_observations*noise_cube_thermal + noise_cube_dets**2)
	#
	logging.info("Saving output")
	if debug_plots:
		plt.clf()

		colors = [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf', u'#2dd42d', u'#eaff00', u'#202020', u'#66f2ff']*2

		fig = plt.figure()
		ax = fig.add_subplot(111)
		
		total_telescope_sky_em = np.zeros_like(lambs_extended)
		
		w, e = np.loadtxt(base_filename + "_sky_em.txt", unpack=True)
		plt.plot(w, e, label="sky", color=colors[-1])
		total_telescope_sky_em += e
		w, e = np.loadtxt(base_filename + "_tel_em.txt", unpack=True)
		plt.plot(w, e, label="telescope", color=colors[-2])
		total_telescope_sky_em += e

		if moon > 0.:
			w, e = np.loadtxt(base_filename + "_moon_em.txt", unpack=True)
			plt.plot(w, e, label="Moon", color=colors[6])
			total_telescope_sky_em += e
		
		
		# HARMONI parts
		total_instrument_em = np.zeros_like(lambs_extended)
		harmoni_files_em = sorted(glob.glob(base_filename + "_HARMONI_*_em.txt"))
		for harmoni_file, color in zip(harmoni_files_em, colors):
			# Read part emission
			w, e = np.loadtxt(harmoni_file, unpack=True)
			m = re.search('.+HARMONI_(.+)_em.txt', harmoni_file)
			plt.plot(w, e, label=m.group(1), color=color, ls="--", lw=1.2)
			# and throuhgput
			w, t = np.loadtxt(base_filename + "_HARMONI_" + m.group(1) + "_tr.txt", unpack=True)
			# Compute the total contribution
			total_instrument_em = total_instrument_em*t + e
		
		plt.plot(w, total_instrument_em, label="HARMONI total", color="red")
		
		np.savetxt(base_filename + "_total_HARMONI_em.txt", np.c_[w, total_instrument_em], comments="#", header="\n".join([
			'TYPE: Total emission.',
			'Wavelength [um], emission']))

		plt.legend(prop={'size': 6})
		plt.xlabel(r"wavelength [$\mu$m]")
		plt.ylabel(r"back emission [photon/m$^2$/$\mu$m/arcsec$^2$]")
		plt.yscale("log")
		#plt.text(0.1, 0.2, "HARMONI/(Telescope+Sky) = {:.2f}".format(np.nanmedian(total_instrument_em/total_telescope_sky_em)), transform=ax.transAxes)
		plt.savefig(base_filename + "_total_em.pdf")

		plt.clf()
		w, e = np.loadtxt(base_filename + "_sky_tr.txt", unpack=True)
		plt.plot(w, e, label="sky", color=colors[-1])
		total_tr_w = w
		total_tr = e

		w, e = np.loadtxt(base_filename + "_tel_tr.txt", unpack=True)
		plt.plot(w, e, label="telescope", color=colors[-2])
		if np.sum(np.abs(total_tr_w - w)) != 0.:
			logging.error('Telescope transmission wavelength error. This should never happen.')
			return
		total_tr *= e
		
		total_instrument_tr = np.ones_like(total_tr)
		# HARMONI parts
		harmoni_files_tr = sorted(glob.glob(base_filename + "_HARMONI_*tr.txt"))
		for harmoni_file, color in zip(harmoni_files_tr, colors):
			w, e = np.loadtxt(harmoni_file, unpack=True)
			m = re.search('.+HARMONI_(.+)_tr.txt', harmoni_file)
			plt.plot(w, e, label=m.group(1), color=color, ls="--", lw=1.2)
			total_instrument_tr *= e
			total_tr *= e
		
		plt.plot(w, total_instrument_tr, label="HARMONI total", color="red")
		
		np.savetxt(base_filename + "_total_HARMONI_tr.txt", np.c_[w, total_instrument_tr], comments="#", header="\n".join([
			'TYPE: Total transmission.',
			'Wavelength [um], transmission']))
		
		# the detector curve has a different wavelength range and spacing
		w, e = np.loadtxt(base_filename + "_det_qe_tr.txt", unpack=True)
		plt.plot(w, e, label="detector", color=colors[5])
		
		total_trans_interp = interp1d(total_tr_w, total_tr, kind='linear', bounds_error=False, fill_value=0.)
		total_tr_final_w = w
		total_tr_final = total_trans_interp(total_tr_final_w)*e

		plt.plot(total_tr_final_w, total_tr_final, label="Total", color=colors[7])
		
		plt.legend(prop={'size': 6})
		plt.xlabel(r"wavelength [$\mu$m]")
		plt.ylabel(r"transmission")
		plt.savefig(base_filename + "_total_tr.pdf")

		if not debug:
			list_files = ["sky_tr", "sky_em", "moon_em", "tel_tr", "tel_em", "det_qe_tr"]
			for _ in list_files:
				os.remove(base_filename + "_" + _ + ".txt")
				os.remove(base_filename + "_" + _ + ".pdf")
				
			for _ in harmoni_files_em+harmoni_files_tr:
				filename = _.replace(".txt", "")
				os.remove(filename + ".txt")
				os.remove(filename + ".pdf")
		
		
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
	outFile_detSNR = base_filename + "_detector_SNR.fits"
	# standard deviation cube
	outFile_std = base_filename + "_std.fits"
	outFile_detstd = base_filename + "_det_std.fits"
	# detectors
	outFile_alldets = base_filename + "_all_dets.fits"
	outFile_useddets = base_filename + "_used_dets.fits"
	
	#
	outFile_read_noise = base_filename + "_read_noise.fits"
	outFile_dark = base_filename + "_dark.fits"
	outFile_ddetector_thermal = base_filename + "_det_thermal.fits"
	
	# Update header
	for _ in simulation_conf:
		head[_.header] = str(_.value)
	
	head['HSM_TIME'] = str(datetime.datetime.utcnow())
	
	# Apply crosstalk to the noiseless cubes
	if grating != "V+R":
		output_cube_spec = apply_crosstalk(output_cube_spec, config_data["crosstalk"])
		output_back_emission_cube = apply_crosstalk(output_back_emission_cube, config_data["crosstalk"])
		output_cube_spec_wo_back = apply_crosstalk(output_cube_spec_wo_back, config_data["crosstalk"])
	
	save_fits_cube(outFile_noiseless_object_plus_back, output_cube_spec*NDIT, "Noiseless O+B", head)
	save_fits_cube(outFile_noiseless_background, output_back_emission_cube*NDIT, "Noiseless B", head)
	save_fits_cube(outFile_noiseless_object, output_cube_spec_wo_back*NDIT, "Noiseless O", head)
	
	save_fits_cube(outFile_observed, sim_total, "Observed O+B1+Noise1", head)
	save_fits_cube(outFile_observed_back, sim_total_only_back, "Observed B2+Noise2", head)
	save_fits_cube(outFile_reduced, sim_reduced, "Reduced (O+B1+Noise1) - (B2+Noise2)", head)
	
	save_fits_cube(outFile_SNR, output_cube_spec_wo_back*NDIT/noise_cube_total, "SNR (O-B)/Noise", head)
	save_fits_cube(outFile_std, noise_cube_total, "Noise standard deviation", head)


	# Flux calibration - from electrons to erg/s/cm2/um
	
	# - Calibration star
	flux_cal_star = 1e-13 # erg/s/cm2/um
	en2ph_conv_fac = (1.98644582e-25 * 1.e7)/(output_lambs*1.e-6*1.e4)
	flux_cal_star_photons = flux_cal_star/en2ph_conv_fac #photon/s/cm2/um

	## r=0.5" aperture for the PSF
	#aperture = CircularAperture((spax_scale.psfsize//2, spax_scale.psfsize//2), r=500./spax_scale.psfscale)
	#psf_fraction = aperture_photometry(psf_internal, aperture)

	flux_cal_star_electrons = flux_cal_star_photons*channel_width*config_data["telescope"]["area"]*np.median(output_transmission) # electron/s  #*psf_fraction['aperture_sum'].data[0] # electron/s
	factor_calibration = flux_cal_star/np.median(flux_cal_star_electrons) # erg/s/cm2/um / (electron/s)
	
	outFile_flux_cal_noiseless = base_filename + "_noiseless_obj_flux_cal.fits"
	outFile_flux_cal_reduced = base_filename + "_reduced_flux_cal.fits"
	
	
	head['BUNIT'] = "erg/s/cm2/um/arcsec2"
	save_fits_cube(outFile_flux_cal_noiseless, output_cube_spec_wo_back/DIT*factor_calibration/spaxel_area, "Flux cal Noiseless O", head)
	save_fits_cube(outFile_flux_cal_reduced, sim_reduced/(NDIT*DIT)*factor_calibration/spaxel_area, "Flux cal Reduced (O+B1+Noise1) - (B2+Noise2)", head)
	
	if det_switch == "True":
		save_fits_cube(outFile_alldets, sim_det_systematics1, "All simulated detectors", head)
		save_fits_cube(outFile_useddets, sim_only_dets, "Used detector noise", head)
		save_fits_cube(outFile_detSNR, output_cube_spec_wo_back*NDIT/noise_cube_total_with_dets, "SNR (O-B)/Exact Noise", head)
		save_fits_cube(outFile_detstd, noise_cube_total_with_dets, "Noise std with exact dets", head)

	if debug:
		save_fits_cube(outFile_dark, noise_cube_dark, "dark noise variance", head)
		save_fits_cube(outFile_read_noise, noise_cube_read_noise**2, "read noise variance", head)
		save_fits_cube(outFile_ddetector_thermal, noise_cube_thermal, "detector thermal noise variance", head)
	
	# Calculate 5-sigma sensitivity
	sens_5sigma = 5.*np.median(noise_cube_total)/DIT*factor_calibration # erg/s/cm2/um
	lcentral = np.median(output_lambs) # micron
	fnu = sens_5sigma*lcentral**2/(const.c.value*1e6) # erg/s/cm2/Hz
	sens_ABmag = -2.5*np.log10(fnu/3631./1e-23) # AB mag
	
	logging.info("Sensitivity point-source 5sigma = {:.2f} mag = {:.2e} erg/s/cm2/um at {:.3f} um".format(sens_ABmag, sens_5sigma, lcentral))

	noise_total = np.median(noise_cube_object + n_observations*noise_cube_back + n_observations*noise_cube_dark + n_observations*noise_cube_thermal + n_observations*noise_cube_read_noise**2)
	fraction_noise_back = n_observations*np.median(noise_cube_back)/noise_total*100.
	fraction_noise_dark = n_observations*np.median(noise_cube_dark)/noise_total*100.
	fraction_noise_read = n_observations*np.median(noise_cube_read_noise)**2/noise_total*100.
	
	logging.info("Noise contributions: Background (sky+tel+instrument) = {:.2f} %. Dark current = {:.2f} %. Read noise = {:.2f} %".format(fraction_noise_back, fraction_noise_dark, fraction_noise_read))


	# Save transmission with crosstalk
	if grating != "V+R":
		output_transmission = apply_crosstalk_1d(output_transmission, config_data["crosstalk"])
	
	np.savetxt(base_filename + "_total_tr.txt", np.c_[output_lambs, output_transmission], comments="#", header="\n".join([
		'TYPE: Total transmission.',
		'Wavelength [um], Transmission']))
	
	
	# Save PSF
	head_PSF = head.copy()
	head_PSF['CDELT1'] = spax_scale.psfscale
	head_PSF['CDELT2'] = spax_scale.psfscale
	head_PSF['LAMBDA'] = psf_lambda
	try:
		del head_PSF['CDELT3']
		del head_PSF['CRVAL3']
		del head_PSF['CRPIX3']
		del head_PSF['CTYPE3']
		del head_PSF['CUNIT3']
		del head_PSF['BUNIT']
		del head_PSF['SPECRES']
	except:
		pass
	
	save_fits_cube(base_filename + "_PSF_internal.fits", psf_internal, "PSF", head_PSF)
	
	# - Rebin and Reshape PSF to match the output cube size and spaxel size
	def rebin_psf(arr, new_shape):
		shape = (new_shape[0], arr.shape[0] // new_shape[0],
			new_shape[1], arr.shape[1] // new_shape[1])
		return arr.reshape(shape).sum(-1).sum(1)
	
	psf_info = config_data["spaxel_scale"][spax]
	
	# Calcualte the offset needed to keep the PSF center
	# at the output image center after rebining
	psf_oversampling = int(round(psf_info.xscale/psf_info.psfscale))
	psf_spaxel_shape = psf_info.psfsize//psf_oversampling + 1
	
	tmp = np.zeros((psf_spaxel_shape*psf_oversampling, psf_spaxel_shape*psf_oversampling))
	psfcenter = psf_info.psfsize//2 + 1
	psfcenter_offset = psfcenter % psf_oversampling
	
	x0 = psf_oversampling//2 + 1 - psfcenter_offset
	tmp[x0:x0 + psf_info.psfsize, x0:x0 + psf_info.psfsize] = psf_internal[:, :]
	
	psf_spaxel = rebin_psf(tmp, (psf_spaxel_shape, psf_spaxel_shape))
	
	if spax == "30x60":
		# an extre rebin is needed for the y axis
		if psf_spaxel_shape % 2 == 1:
			tmp = np.zeros((psf_spaxel_shape + 1, psf_spaxel_shape + 1))
			final_psf_shape = ((psf_spaxel_shape + 1)//2, psf_spaxel_shape + 1)
			tmp[1:, 1:] = psf_spaxel
		else:
			tmp = psf_spaxel
			final_psf_shape = (psf_spaxel_shape//2, psf_spaxel_shape)
		
		psf_spaxel = rebin_psf(tmp, final_psf_shape)
	
	# center PSF on the output array
	center_x_output = (output_cube_spec.shape[2] - 1) // 2 - (output_cube_spec.shape[2] % 2 - 1)
	center_y_output = (output_cube_spec.shape[1] - 1) // 2 - (output_cube_spec.shape[1] % 2 - 1)
	center_x_psf = (psf_spaxel.shape[1] - 1) // 2 - (psf_spaxel.shape[1] % 2 - 1)
	center_y_psf = (psf_spaxel.shape[0] - 1) // 2 - (psf_spaxel.shape[0] % 2 - 1)

	# adjust y axis
	tmp = np.zeros((output_cube_spec.shape[1], psf_spaxel.shape[1]))
	if tmp.shape[0] > psf_spaxel.shape[0]:
		yi = center_y_output - center_y_psf
		tmp[yi:yi+psf_spaxel.shape[0], :] = psf_spaxel
	else:
		yi = center_y_psf - center_y_output
		tmp[:, :] = psf_spaxel[yi:yi+tmp.shape[0], :]
	psf_spaxel = tmp

	# adjust x axis
	tmp = np.zeros((output_cube_spec.shape[1], output_cube_spec.shape[2]))
	if tmp.shape[1] > psf_spaxel.shape[1]:
		xi = center_x_output - center_x_psf
		tmp[:, xi:xi+psf_spaxel.shape[1]] = psf_spaxel
	else:
		xi = center_x_psf - center_x_output
		tmp[:, :] = psf_spaxel[:, xi:xi+tmp.shape[1]]
	psf_spaxel = tmp
	
	head_PSF['CDELT1'] = spax_scale.xscale
	head_PSF['CDELT2'] = spax_scale.yscale
	save_fits_cube(base_filename + "_PSF.fits", psf_spaxel, "PSF", head_PSF)
	
	


	if hsimlog.count_error == 0 and hsimlog.count_warning == 0:
		logging.info('Simulation OK - ' + str(hsimlog.count_error) + " errors and " + str(hsimlog.count_warning) + " warnings")
	else:
		logging.warning('Simulation with problems - ' + str(hsimlog.count_error) + " errors and " + str(hsimlog.count_warning) + " warnings")
	
	
	return
	

def save_fits_cube(filename, data, typ, header):
	header['HSM_TYPE'] = typ
	fits.writeto(filename, data, header=header, overwrite=True, output_verify="silentfix")
	
	

class HSIMLoggingHandler(logging.Handler):
	def __init__(self):
		self.count_warning = 0
		self.count_error = 0
		logging.Handler.__init__(self)
		
	def emit(self, record):
		if record.levelno == logging.WARNING:
			self.count_warning += 1
		elif record.levelno >= logging.ERROR:
			self.count_error += 1
		

class HSIMFormatter(logging.Formatter):

	def __init__(self):
		logging.Formatter.__init__(self, "%(message)s")

	def format(self, record):
		if record.levelno > logging.INFO:
			self._fmt = "** %(levelname)s ** %(message)s"
		else:
			self._fmt = "%(message)s"

		return logging.Formatter.format(self, record)
