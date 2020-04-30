'''
Misc utilities
'''

import os
import logging

#Datafile path setup
def path_setup(path):
	'''function to setup paths to datafiles'''
	script_dir = os.path.split(os.path.abspath(__file__))[0]
	out_path = os.path.join(script_dir, os.path.relpath(path),'')

	return out_path

def trim_cube(datacube, verbose=False):
	''' Trims datacube to a size that will fit on 8 HAWAII 4RG
		detectors when detector systematics are used.

	Inputs:
		datacube: object datacube

	Outputs:
		new_cube: trimmed cube

	'''
	z, x, y = datacube.shape
	if z > 3700:
		lam_low = int((z - 3700) / 2)
		if z % 2 == 0:
			lam_high = int(z - lam_low)
		else:
			lam_high = int(z - lam_low - 1)
		if verbose:
			logging.info('Restricting spectral axis to 3700 channels')
	else:
		lam_low, lam_high = 0, z

	if x > 152:
		x_low = int((x - 152) / 2)
		if x % 2 == 0:
			x_high = int(x - x_low)
		else:
			x_high = int(x - x_low - 1)
		if verbose:
			logging.info('Restricting x axis to 152 channels')
	else:
		x_low, x_high = 0, x
		
	if y > 204:
		y_low = int((y - 204) / 2)
		if y % 2 == 0:
			y_high = int(y - y_low)
		else:
			y_high = int(y - y_low - 1)
		if verbose:
			logging.info('Restricting y axis to 204 spaxels')
	else:
		y_low, y_high = 0, y
		
	new_cube = datacube[lam_low:lam_high, x_low:x_high, y_low:y_high]
	return new_cube
