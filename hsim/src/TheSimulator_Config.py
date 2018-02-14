'''File that stores hardwired data for use in HARMONI
simulation pipeline. Data stored in dictionary format
with keywords.

Written by Simon Zieleniewski

Last updated 27-04-16

'''

config_data = {

    'E-ELT': {'diam':37., 'obsc':0.3, 'area':932.46},#diam [m], area [m^2]
#	'E-ELT': {'diam':37., 'obsc':0.3, 'area':875.6},#diam [m], area [m^2]

    'VLT': {'diam':8.2, 'obsc':0.1146, 'area':52.1},

    'emissivity':0.244,

    'dark_current_vis':0.00042,     #e/pix/s

    'dark_current_nir':0.0053,      #e/pix/s

    'read_noise_vis':2.0,           #e/pix

    'read_noise_nir':2.845,         #e/pix

    'read_noise_nir_lowexp':12.0,   #e/pix

    'trans_wo_grat':0.35,     #Default=0.35, HARMONI transmission without grating curve

    'ADR': 'ON',                     #ADR effect ON or OFF

    #HARMONI PDR gratings (lambda_min, lambda_max, R)
    'gratings': {'V+R':(0.4625,0.8200, 3100.), 
		'Iz+J':(0.8105,1.3695, 3300.), 'H+K':(1.4500,2.4500, 3300.),
		'Iz':(0.8300,1.0500, 7100.),'J':(1.0463,1.3237, 7100.), 'H':(1.4348,1.8151, 7100.), 'K':(1.9514,2.4686, 7100.),
		'z':(0.8268,0.9032, 17200.), 'J-high':(1.1900,1.3000, 17200.),
		'H-high':(1.5341,1.6759, 17200.), 'K-high':(2.0933,2.2867, 17200.),
		'None':None},

    'field_of_view': (152, 214),    #(x, y) field of view in spaxels


    }
