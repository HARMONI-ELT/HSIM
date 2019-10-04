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

    'trans_w_grat':0.45,     #Default=0.45, HARMONI transmission when using grating curve
    'trans_wo_grat':0.35,     #Default=0.35, HARMONI transmission without grating curve

    'ADR': 'ON',                     #ADR effect ON or OFF

    #HARMONI phase A gratings (lambda_min, lambda_max, R)
    'gratings': {'V':(.47,.63, 7500.),'R':(.63,.79, 7500.), 'Iz':(.82,1.03, 7500.),
             'J':(1.08,1.36, 7500.), 'H':(1.46,1.83, 7500.), 'K':(1.95,2.45, 7500.),
             'V+R':(.47,.81, 3500.), 'Iz+J':(.8,1.36, 3500.), 'H+K':(1.45,2.45, 3500.),
             'V-high':(.53,.59, 20000.), 'R-high':(.61,.68, 20000.),
             'z':(.82,.91, 20000.), 'J-high':(1.17,1.29, 20000.),
             'H-high':(1.545,1.715, 20000.), 'K-high':(2.09,2.32, 20000.),
             'None':None},

    'field_of_view': (152, 214),    #(x, y) field of view in spaxels


    }
