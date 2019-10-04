'''File that stores hardwired data for use in JWST/NIRSpec
simulation pipeline. Data stored in dictionary format
with keywords.

Written by Simon Zieleniewski

Created: 16-10-15

Last updated 03-12-15

'''

config_data = {

    'area':25.911,                  #m^2

    'emissivity':0.15,              #dimensionless (not used)

    'dark_current':0.01,            #e/pix/s

    'read_noise':5.1,               #e/pix

    'transmission':'intrinsic',     #Default=intrinsic, or enter value here if desired

    #NIRSpec gratings (lambda_min, lambda_max, delta_lambda)
    'gratings': {'F070-G140M':(.7,1.2,0.00140), 'F100-G140M':(1.0,1.8,0.00140),
                 'G235M':(1.7,3.0,0.00237), 'G395M':(2.9,5.0,0.00400),
                 'F070-G140H':(0.7,1.2,0.000516), 'F100-G140H':(1.0,1.8,0.000516),
                 'G235H':(1.7,3.0,0.00087), 'G395H':(2.9,5.0,0.000146),
                 'Prism':'Lowres'},

    'spaxels': (100.,100.),         #Spaxel scale [mas, mas]

    'field_of_view': (30, 30),      #(x, y) field of view in spaxels


    }
