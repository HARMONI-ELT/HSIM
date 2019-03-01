'''
Code for the bulk framework of the HARMONI
simulator. This code should show the main functions/processes
to move from an input datacube (lambda, y, x) to output cubes:
'''
import collections
import datetime
import multiprocessing as mp
import os.path
import logging

import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import convolve2d
from astropy.io import fits

import matplotlib.pylab as plt

from src.config import *
from src.init_cube import init_cube
from src.sim_sky import sim_sky
from src.sim_telescope import sim_telescope
from src.sim_instrument import sim_instrument
from src.sim_detector import sim_detector, apply_crosstalk, mask_saturated_pixels, apply_crosstalk_1d
from src.sim_detector import make_det_instance, add_detectors
from src.modules.adr import apply_adr
from src.modules.rebin import *
from src.modules.misc_utils import trim_cube


def main(datacube, outdir, DIT, NDIT, grating, spax, seeing, air_mass, version, res_jitter=3., moon=0.,
	 site_temp=280.5, adr_switch='True', det_switch='False', seednum=100, nprocs=mp.cpu_count()-1, debug=False, aoMode="LTAO"):
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
		det_switch: Boolean - use detector systematics (off by default)
		seednum: ramdom seed number
		nprocs: Number of processes
		debug: keep debug plots and all noise outputs
		aoMode: Adaptive optics mode: "LTAO", "SCAO", or "noAO" for seeing limited

	Outputs:

	'''
	debug_plots = True
	aoMode = aoMode.upper()

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
			Conf('AO', 'HSM_AO', aoMode),
			Conf('No. of processes', 'HSM_NPRC', nprocs),
			]

	base_name = os.path.splitext(os.path.basename(datacube))[0]
	base_filename = os.path.join(outdir, base_name)

	logfile = base_filename + ".log"
	open(logfile, 'w').close()

	# log to a file
	logging.basicConfig(filename=logfile, level=logging.DEBUG, format='%(asctime)s  %(levelname)s  %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
	# and also to the terminal
	logger = logging.getLogger()
	std = logging.StreamHandler()
	std.setFormatter(HSIMFormatter())
	logger.addHandler(std)

	hsimlog = HSIMLoggingHandler()
	logger.addHandler(hsimlog)


	logging.info("Simulation input parameters:")
	for _ in simulation_conf:
		logging.info(_.name + " = " + str(_.value))

	if aoMode not in ["LTAO", "SCAO", "NOAO", "AIRY"]:
		logging.error(aoMode + ' is not a valid AO mode. Valid options are: LTAO, SCAO, noAO, Airy')
		return

	if air_mass not in config_data["PSD_cube"]["air_masses"]:
		logging.error(str(air_mass) + ' is not a valid air mass. Valid options are: ' + ", ".join(map(str, sorted(config_data["PSD_cube"]["air_masses"]))))
		return

	if seeing not in  config_data["PSD_cube"]["seeings"]:
		logging.error(str(seeing) + ' is not a valid seeing. Valid options are: ' + ", ".join(map(str, sorted(config_data["PSD_cube"]["seeings"]))))
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
	cube_exp, back_emission, transmission, psf, psf_lambda = sim_telescope(cube_exp, back_emission, transmission, lambs_extended, cube_lamb_mask, DIT, res_jitter, air_mass, seeing, spax, site_temp, aoMode, nprocs, debug_plots=debug_plots, output_file=base_filename)

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
	
	output_cube_spec = np.zeros((len(output_lambs), out_size_y, out_size_x))
	
	logging.info("- Output spectral sampling: {:.2f} A".format(new_lamb_per_pix*10000.))
	
	for i in np.arange(0, out_size_x):
		for j in np.arange(0, out_size_y):
			output_cube_spec[:, j, i] = rebin1d(output_lambs, lambs, output_cube[:, j, i])
	
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
	spaxel_area = spax_scale.xscale/1000.*spax_scale.yscale/1000. # [arcsec2]
	channel_width = new_lamb_per_pix
	
	output_cube_spec = output_cube_spec*spaxel_area*channel_width*config_data["telescope"]["area"]
	output_back_emission = output_back_emission*spaxel_area*channel_width*config_data["telescope"]["area"]
	head['FUNITS'] = "photons"
	
	# 6 - Detector
	#	- QE
	#	- Dark
	#	- Read noise

	# Cut cubes to correct size if using detector systematics and generate detectors 
	if det_switch == "True":
		logging.info("Trimming datacubes to correct size")
		print output_cube_spec.shape
		output_cube_spec = trim_cube(output_cube_spec)
		print output_cube_spec.shape
		logging.info("Generating simulated detectors")
		sim_dets1 = make_det_instance(NDIT)
		sim_dets2 = make_det_instance(NDIT)
	
	output_cube_spec, output_back_emission, output_transmission, read_noise, dark_current = sim_detector(output_cube_spec, output_back_emission, output_transmission, output_lambs, grating, DIT, debug_plots=debug_plots, output_file=base_filename)

	head['FUNITS'] = "electrons"
	
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
		
	if det_switch == "False":
		# - read noise and dark current for object exposure
		sim_read_noise1 = np.random.normal(zero_cube, np.sqrt(NDIT)*read_noise).astype(np.float32)
		sim_dark_current1 = np.random.poisson(dark_cube).astype(np.float32)
	
		# - combine object, read noise and dark current
		sim_total = sim_object_plus_back + sim_read_noise1 + sim_dark_current1
	else:
		sim_total = add_detectors(sim_object_plus_back, sim_dets1)
	
	# B. Background exposure
	#- background with crosstalk
	sim_back = np.random.poisson(abs(output_back_emission_cube*NDIT)).astype(np.float32)
	# Apply crosstalk only to NIR detectors
	if grating != "V+R":
		sim_back = apply_crosstalk(sim_back, config_data["crosstalk"])


	if np.sum(saturated_back) > 0:
		logging.warning(str(np.sum(saturated_back)) + " pixels are saturated in the back frames")
		sim_back[saturated_back] = np.nan

	if det_switch == "False":
		# - read noise and dark current for background exposure
		sim_read_noise2 = np.random.normal(zero_cube, np.sqrt(NDIT)*read_noise).astype(np.float32)
		sim_dark_current2 = np.random.poisson(dark_cube).astype(np.float32)
	
		# - combine object, read noise and dark current
		sim_total_only_back = sim_back + sim_read_noise2 + sim_dark_current2
	else:
		sim_total_only_back = add_detectors(sim_back, sim_dets2)
	
	
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
	
	noise_cube_total = np.sqrt(noise_cube_object + 2.*noise_cube_back + 2.*noise_cube_dark + 2.*noise_cube_read_noise**2)
	
	#
	logging.info("Saving output")
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
		
		if aoMode not in ["NOAO", "AIRY"]:
			w, e = np.loadtxt(base_filename + "_ins_AOd_em.txt", unpack=True)
			plt.plot(w, e, label="AO dichroic", color=colors[2])
			total_instrument_em += e
		
		w, e = np.loadtxt(base_filename + "_ins_FPRS_em.txt", unpack=True)
		plt.plot(w, e, label="FPRS", color=colors[3])
		total_instrument_em += e
		w, e = np.loadtxt(base_filename + "_ins_em.txt", unpack=True)
		plt.plot(w, e, label="instrument after FPRS", color=colors[4])
		total_instrument_em += e
		
		plt.legend(prop={'size': 6})
		plt.xlabel(r"wavelength [$\mu$m]")
		plt.ylabel(r"back emission [photons/m$^2$/$\mu$m/arcsec$^2$]")
		plt.yscale("log")
		plt.text(0.1, 0.2, "HARMONI/(Telescope+Sky) = {:.2f}".format(np.nanmedian(total_instrument_em/total_telescope_sky_em)), transform=ax.transAxes)
		plt.savefig(base_filename + "_total_em.pdf")

		plt.clf()
		w, e = np.loadtxt(base_filename + "_sky_tr.txt", unpack=True)
		plt.plot(w, e, label="sky", color=colors[0])
		total_tr_w = w
		total_tr = e

		w, e = np.loadtxt(base_filename + "_tel_tr.txt", unpack=True)
		plt.plot(w, e, label="telescope", color=colors[1])
		if np.sum(np.abs(total_tr_w - w)) != 0.:
			logging.error('Telescope transmission wavelength error. This should never happen.')
			return
		total_tr *= e
	
		if aoMode not in ["NOAO", "AIRY"]:
			w, e = np.loadtxt(base_filename + "_ins_AOd_tr.txt", unpack=True)
			plt.plot(w, e, label="AO dichroic", color=colors[2])
			if np.sum(np.abs(total_tr_w - w)) != 0.:
				logging.error('AO transmission wavelength error. This should never happen.')
				return
			total_tr *= e


		w, e = np.loadtxt(base_filename + "_ins_FPRS_tr.txt", unpack=True)
		plt.plot(w, e, label="FPRS", color=colors[3])
		if np.sum(np.abs(total_tr_w - w)) != 0.:
			logging.error('FPRS transmission wavelength error. This should never happen.')
			return
		total_tr *= e
		
		w, e = np.loadtxt(base_filename + "_ins_tr.txt", unpack=True)
		plt.plot(w, e, label="instrument after FPRS", color=colors[4])
		if np.sum(np.abs(total_tr_w - w)) != 0.:
			logging.error('instrument transmission wavelength error. This should never happen.')
			return
		total_tr *= e
		
		# the detector curve has a different wavelength range and spacing
		w, e = np.loadtxt(base_filename + "_det_qe.txt", unpack=True)
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
			list_files = ["sky_tr", "sky_em", "moon_em", "tel_tr", "tel_em", "ins_tr", "ins_em", "det_qe", "ins_FPRS_tr", "ins_FPRS_em"]
			if aoMode not in ["NOAO", "AIRY"]:
				list_files.append("ins_AOd_tr")
				list_files.append("ins_AOd_em")
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
	# standard deviation cube
	outFile_std = base_filename + "_std.fits"
	
	# 
	outFile_read_noise = base_filename + "_read_noise.fits"
	outFile_dark = base_filename + "_dark.fits"
	
	# Update header
	for _ in simulation_conf:
		head[_.header] = str(_.value)
	
	head['HSM_TIME'] = str(datetime.datetime.utcnow())
	
	# Apply to the noiseless cubes
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
	
	if debug:
		save_fits_cube(outFile_dark, noise_cube_dark, "dark noise variance", head)
		save_fits_cube(outFile_read_noise, noise_cube_read_noise**2, "read noise variance", head)
	

	# Save transmission with crosstalk
	if grating != "V+R":
		output_transmission = apply_crosstalk_1d(output_transmission, config_data["crosstalk"])
	
	np.savetxt(base_filename + "_total_tr.txt", np.c_[output_lambs, output_transmission], comments="#", header="\n".join([
		'TYPE: Total transmission.',
		'Wavelength [um], Transmission']))
	
	
	# Save PSF
	head_PSF = head
	head_PSF['CDELT1'] = spax_scale.psfscale
	head_PSF['CDELT2'] = spax_scale.psfscale
	head_PSF['LAMBDA'] = psf_lambda
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
	
	if hsimlog.count_error == 0 and hsimlog.count_warning == 0:
		logging.info('Simulation OK')
	else:
		logging.warning('Simulation with problems - ' + str(hsimlog.count_error) + " errors and " + str(hsimlog.count_warning) + " warnings")
	
	
	return
	

def save_fits_cube(filename, data, typ, header):
	header['HSM_TYPE'] = typ
	fits.writeto(filename, data, header=header, overwrite=True)
	
	

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
