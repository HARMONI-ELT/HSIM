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
	'detector':{
		"avg":{
			'read_noise': {"vis":2.0, "nir":6.0, "nir_lowexp":15.0}, # e/pix
			'dark_current': {"vis":0.00042, "nir":0.0053}, # e/pix/s
			'qe':'file',
		},
		"best":{
			'read_noise': {"vis":2.0, "nir":4.0, "nir_lowexp":10.0}, # e/pix
			'dark_current': {"vis":0.00042, "nir":0.0042}, # e/pix/s
			'qe':{"w":[0.9, 1.2, 2.0], "qe":[0.85, 0.90, 0.95]},
		},
		"worst":{
			'read_noise': {"vis":2.0, "nir":8.0, "nir_lowexp":12.0}, # e/pix
			'dark_current': {"vis":0.00042, "nir":0.0220}, # e/pix/s
			'qe':{"w":[0.9, 1.2, 2.0], "qe":[0.70, 0.80, 0.90]},
		},
		"contractual":{
			'read_noise': {"vis":2.0, "nir":30.0, "nir_lowexp":30.0}, # e/pix
			'dark_current': {"vis":0.00042, "nir":0.1}, # e/pix/s
			'qe':{"w":[0.9, 1.2, 2.0], "qe":[0.50, 0.50, 0.50]},
		},
	},
	
	'saturation': {"vis":72000., "nir":30000.}, # e
	
	'crosstalk': 0.02, # fraction of photons going to each of the 4 contiguous pixels
	'side_length':4096,
	'N_IR_det':8,

	#Detector systematics parameters
	'systematics': {"rd":2.845,
                        "rd_lowexp":12.0,
                        "rn_file":"kmos_rn.fits",
                        "pedestal":4,
                        "c_pink":3,
                        "u_pink":1,
                        "acn":0.5,
                        "pca0_amp":0.2,
                        "ref_ratio":0.8,
                        "ktc_noise":29.,
                        "bias_offset":5000.,
                        "bias_amp":500.,
                        "force_new":False
                        },
                        
	'spectral_sampling':{"output":2.2, "internal":4.}, # spectral sampling of the output cube and internal. Nyquist = 2
	'LSF_kernel_size':12., # LSF kernel size in sigma units
	
	'telescope': {'diameter':37., 'area':{"VRIz":996.3, "JHK":896.3}}, #diam [m], area [m^2]
	
	'HARMONI_FPRS_temp': +2., # C
	'HARMONI_cryo_temp': 130., # K
	
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
	
	'gratings_nominal': {	#low resolution
			'V+R':GratingInfo(0.458, 0.8200, 3000.), ##
			'Iz+J':GratingInfo(0.83, 1.09, 3000.),
			'H+K':GratingInfo(1.450, 1.925, 3000.),
			# med-resolution
			'Iz':GratingInfo(0.830, 0.94, 7000.),
			'J':GratingInfo(1.05, 1.185, 7000.),
			'H':GratingInfo(1.45, 1.625, 7000.),
			'K':GratingInfo(1.97, 2.185, 7000.),
			# high-resolution
			'z-high':GratingInfo(0.828, 0.865, 17000.),
			'J-short':GratingInfo(1.012, 1.102, 17000.), ##
			'J-long':GratingInfo(1.098, 1.189, 17000.), ##
			'H-high':GratingInfo(1.538, 1.608, 17000.),
			'K-short':GratingInfo(2.017, 2.201, 17000.),
			'K-long':GratingInfo(2.199, 2.300, 17000.)
			},
	
	'spaxel_scale': {'4x4':SpaxelScaleInfo(4., 4., 0.8, 1250),
		  '7x7':SpaxelScaleInfo(7., 7., 1.4, 1250),
		  '10x10':SpaxelScaleInfo(10., 10., 2., 800),
		  '20x20':SpaxelScaleInfo(20., 20., 4., 580),
		  '25x25':SpaxelScaleInfo(25., 25., 5., 580),
		  '30x60':SpaxelScaleInfo(30., 60., 6., 400),
		  '60x60':SpaxelScaleInfo(60., 60., 6., 400),
		  '120x60':SpaxelScaleInfo(120., 60., 6., 400)
		  },
	
	
	#FWHM of Instrument PSF depending on output spaxel scale in mas
	#Factors taken into account:
	#design image quality, manufacturing and assembly tolerances, vibration, flexure, diff refraction,
	'dynamic_instrument_psf': 5.5,
	'static_instrument_psf': {'4x4': 3.,
		  '7x7': 7.,
    		  '10x10':14.,
		  '20x20':28.,
		  '25x25':12.,
		  '30x60':30.,
		  '60x60':30.,
		  '120x60':30.
		},
	
	# minimum compliant instrument
	'mci_dynamic_instrument_psf': 5.5,
	'mci_static_instrument_psf': {'4x4': 3.,
		  '7x7':7.,
		  '10x10':14.,
		  '20x20':28.,
		  '25x25':12.,
		  '30x60':135.,
		  '60x60':135.,
		  '120x60':135.
		},

	#Each PSD file containts 1 seeing  [0.43] and 1 zenith angle [25]
	'PSD_file':{"LTAO":"psd_ltao_hsim_6LGS_cn2_310.fits", 
		"SCAO":"psd_SCAO_hsim_6_cn2_310.fits"},
	'PSD_params':{'air_masses':[1.1, 1.3, 1.5, 2.0],
		'seeings':[0.43, 0.57, 0.64, 0.72, 1.04]}

}
	
	
def get_telescope_area(grating_name):
	if grating_name in ['V+R', 'Iz+J', 'Iz', 'z-high']:
		return config_data["telescope"]["area"]["VRIz"]
	else:
		return config_data["telescope"]["area"]["JHK"]
	


class HSIMError(Exception):
	pass
	def __init__(self, message):
		logging.error(message)
		sys.exit()

