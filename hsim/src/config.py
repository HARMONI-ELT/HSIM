'''File that stores hardwired data for use in HARMONI
simulation pipeline. Data stored in dictionary format
with keywords.
'''
import logging
import collections

GratingInfo = collections.namedtuple('GratingInfo', 'lmin, lmax, R')
SpaxelScaleInfo = collections.namedtuple('SpaxelScaleInfo', 'xscale, yscale, psfscale, psfsize')
	
config_data = {
	'read_noise': {"vis":2.0, "nir":2.845, "nir_lowexp":12.0}, # e/pix
	'dark_current': {"vis":0.00042, "nir":0.0053}, # e/pix/s
	'saturation': {"vis":55000., "nir":55000.}, # e
	'crosstalk': 0.02, # fraction of photons going to each of the 4 contiguous pixels

	'spectral_sampling':{"output":2.2, "internal":4.}, # spectral sampling of the output cube and internal. Nyquist = 2
	'LSF_kernel_size':12., # LSF kernel size in sigma units
	
	'telescope': {'diameter':37., 'obsc':0.3, 'area':932.46}, #diam [m], area [m^2]
	
	'HARMONI_temp':140., # K
	'HARMONI_FPRS_diff_temp':20., # Tamb - 20 K for the FPRS
	'HARMONI_AO_diff_temp':5., #  Tamb - 5 K for the AO dichroic
	
	'data_dir':"sim_data/",
		
	'gratings': {'V+R':GratingInfo(0.4625,0.8200, 3100.), 
			'Iz+J':GratingInfo(0.8105,1.3695, 3300.), 'H+K':GratingInfo(1.4500,2.4500, 3300.),
			'Iz':GratingInfo(0.8300,1.0500, 7100.),'J':GratingInfo(1.0463,1.3237, 7100.), 'H':GratingInfo(1.4348,1.8151, 7100.), 'K':GratingInfo(1.9514,2.4686, 7100.),
			# high-resolution
			'z-high':GratingInfo(0.829, 0.902, 17200.),
			'J-short':GratingInfo(1.012, 1.102, 17200.),
			'J-long':GratingInfo(1.098, 1.189, 17200.),
			'J-high':GratingInfo(1.1900, 1.3000, 17200.),
			'H-high':GratingInfo(1.541, 1.676, 17200.),
			'K-short':GratingInfo(2.024, 2.202, 17200.),
			'K-long':GratingInfo(2.197, 2.390, 17200.)
			},
	
	'spaxel_scale': {'4x4':SpaxelScaleInfo(4., 4., 1., 1000),
		  '10x10':SpaxelScaleInfo(10., 10., 2., 800),
		  '20x20':SpaxelScaleInfo(20., 20., 3., 600),
		  '30x60':SpaxelScaleInfo(30., 60., 3., 600)
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
		logging.error(message)

