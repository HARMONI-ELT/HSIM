'''File that stores hardwired data for use in HARMONI
simulation pipeline. Data stored in dictionary format
with keywords.
'''
import logging
import collections
import sys

GratingInfo = collections.namedtuple('GratingInfo', 'lmin, lmax, R')
SpaxelScaleInfo = collections.namedtuple('SpaxelScaleInfo', 'xscale, yscale, psfscale, psfsize')
	
config_data = {
	'read_noise': {"vis":2.0, "nir":2.845, "nir_lowexp":12.0}, # e/pix
	'dark_current': {"vis":0.00042, "nir":0.0053}, # e/pix/s
	'saturation': {"vis":55000., "nir":55000.}, # e
	'crosstalk': 0.02, # fraction of photons going to each of the 4 contiguous pixels
	'side_length': 4096, # length of side of H4RG

	'spectral_sampling':{"output":2.2, "internal":4.}, # spectral sampling of the output cube and internal. Nyquist = 2
	'LSF_kernel_size':12., # LSF kernel size in sigma units
	
	'telescope': {'diameter':37., 'area':980.0}, #diam [m], area [m^2]
	
	'HARMONI_FPRS_diff_temp':20., # Tamb - 20 K for the FPRS
	'HARMONI_cryo_temp':130., # K
	
	'data_dir':"sim_data/",
	
	#  HRM-00244
	'gratings': {	#low resolution
			'V+R':GratingInfo(0.458, 0.8200, 3100.),
			'Iz+J':GratingInfo(0.811, 1.369, 3355.),
			'H+K':GratingInfo(1.450, 2.450, 3355.),
			# med-resolution
			'Iz':GratingInfo(0.830, 1.050, 7104.),
			'J':GratingInfo(1.046, 1.324, 7104.),
			'H':GratingInfo(1.435, 1.815, 7104.),
			'K':GratingInfo(1.951, 2.469, 7104.),
			# high-resolution
			'z-high':GratingInfo(0.828, 0.902, 17385.),
			'J-short':GratingInfo(1.012, 1.102, 17385.),
			'J-long':GratingInfo(1.098, 1.189, 17385.),
			'H-high':GratingInfo(1.538, 1.678, 17385.),
			'K-short':GratingInfo(2.017, 2.201, 17385.),
			'K-long':GratingInfo(2.199, 2.400, 17385.)
			},
	
	'spaxel_scale': {'4x4':SpaxelScaleInfo(4., 4., 0.8, 1250),
		  '10x10':SpaxelScaleInfo(10., 10., 2., 800),
		  '20x20':SpaxelScaleInfo(20., 20., 2.857, 800),
		  '30x60':SpaxelScaleInfo(30., 60., 3.333, 600)
		  },
	
	
	#FWHM of Instrument PSF depending on output spaxel scale in mas
	#Factors taken into account:
	#design image quality, manufacturing and assembly tolerances, vibration, flexure, diff refraction,
	'dynamic_instrument_psf': 4.8,
	'static_instrument_psf': {'4x4': 0.,
    		  '10x10':10.,
		  '20x20':26.,
		  '30x60':30.
		},

	#Each PSD file containts 5 seeings [0.43, 0.57, 0.64, 0.72, 1.04] and 4 zenith angles [25, 40, 48, 60]
	'PSD_file':{"LTAO":"psd_ltao_hsim_6LGS_cn2_all.fits", 
		"SCAO":"psd_SCAO_hsim_6_cn2_all.fits"},
	'PSD_cube':{'air_masses':[1.1, 1.3, 1.5, 2.0],
		'seeings':[0.43, 0.57, 0.64, 0.72, 1.04]}

}

class HSIMError(Exception):
	pass
	def __init__(self, message):
		sys.tracebacklimit = 0
		logging.error(message)

