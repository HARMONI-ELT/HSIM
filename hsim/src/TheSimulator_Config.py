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
    'gratings': {'V+R':(.47,.81, 3500.), 'Iz+J':(.83,1.37, 3500.), 'H+K':(1.43,2.45, 3500.),
	    'Iz':(.85,1.00, 7500.),'J':(1.03,1.37, 7500.), 'H':(1.43,1.82, 7500.), 'K':(1.93,2.45, 7500.),
             'z':(.825,.900, 17000.), 'J-high':(1.175,1.325, 17000.),
             'H-high':(1.525,1.675, 17000.), 'K-high':(2.075,2.275, 17000.),
             'None':None},

    'field_of_view': (152, 214),    #(x, y) field of view in spaxels


    }
