'''
Misc utilities
'''

import os

#Datafile path setup
def path_setup(path):
	'''function to setup paths to datafiles'''
	script_dir = os.path.split(os.path.abspath(__file__))[0]
	out_path = os.path.join(script_dir, os.path.relpath(path),'')

	return out_path

