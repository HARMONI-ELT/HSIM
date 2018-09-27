'''
Reads input FITS data cube and resamples according to the
grating and spaxel scale. The output cube is in ph/s/m2/um/arcsec2 units
'''
import os
import logging

import astropy.io.fits as fits
import numpy as np
import scipy.constants as sp
from scipy.interpolate import interp2d

from modules.fits_utils import *
from config import *
from modules.rebin import *

def spectral_res(datacube, head, grating, wavels):
	'''Function that takes input datacube and rebins it to the
	chosen spectral resolution. It interpolates all spaxels along wavelength axis
	then extracts pixel values for each datacube wavelength channel. Function combines both
	spectral resolution choice and datacube chopping in wavelength axis.

	Inputs:
		datacube: object datacube
		head: header file for object datacube
		grating: string representing chosen grating.
		wavels: wavelength array for datacube

	Outputs:
		new_cube: rebinned datacube
		head: updated header file
		new_wavels: new wavelength array for datacube
		input_spec_res: input spectral resolution [micron]
	'''

	logging.info('# Spectral resolution and wavelength range')
	logging.info('Chosen grating: ' + grating)
	z, y, x = datacube.shape
	
	try:
		bandws = config_data['gratings'][grating]
		new_res = (bandws.lmin + bandws.lmax)/(2.*bandws.R)
		lamb_per_pix = new_res/config_data["spectral_sampling"]["internal"]
	except:
		raise HSIMError(grating + ' is not a valid grating. Valid options are: ' + ", ".join(sorted(config_data['gratings'].keys())))
	
	input_spec_res = head['SPECRES'] # micron
	if input_spec_res > new_res:
		logging.warning("The input spectral resolution is lower than the HARMONI grating resolution. Assuming input resolution = HARMONI resolution")
		input_spec_res = new_res
	elif input_spec_res > 0.5*new_res:
		logging.warning("The input spectral resolution is lower than 2 times the HARMONI grating resolution")

	
	logging.info('Input resolution = %.1f A' % (input_spec_res*10000.))
	logging.info('Input sampling = %.1f A' % (head['CDELT3']*10000.))
	logging.info('Output resolution = %.1f A' % (new_res*10000.))
	logging.info('Internal sampling = %.1f A' % (lamb_per_pix*10000.))

	#Interpolate datacube onto regular pixel grid
	if wavels[0] <= bandws.lmin and wavels[-1] >= bandws.lmax: # cube wavelength range larger than grating choice
		logging.warning('Cube wavelength range larger than grating: chopping down to grating range')
		new_wavels = np.arange(bandws.lmin, bandws.lmax, lamb_per_pix)
		
	elif wavels[0] > bandws.lmin and wavels[-1] < bandws.lmax: # cube wavelength range inside grating choice
		logging.info('Cube wavelength range within grating range')
		new_wavels = np.arange(wavels[0], wavels[-1], lamb_per_pix)
	
	elif wavels[0] < bandws.lmin and wavels[-1] < bandws.lmax and wavels[-1] > bandws.lmin: # cube short wavelength longer than grating shortest wavelength
		logging.warning('Cube shortest wavelength shorter than grating')
		new_wavels = np.arange(bandws.lmin, wavels[-1], lamb_per_pix)
		
	elif wavels[0] > bandws.lmin and wavels[0] < bandws.lmax and wavels[-1] > bandws.lmax: # cube short wavelength longer than grating shortest wavelength
		logging.warning('Cube longest wavelength longer than grating')
		new_wavels = np.arange(wavels[0], bandws.lmax, lamb_per_pix)
		
	else:
		raise HSIMError('The wavelength range of the input cube ' + str(wavels[0]) + " - " + str(wavels[1]) + " is not valid for the " + grating + " grating ("  + str(bandws.lmin) + " - " + str(bandws.lmax) + ")")

	new_cube = np.zeros((len(new_wavels), y, x), dtype=float)
	
	# regrid spectrum and conserve flux
	logging.info('Interpolating data cube - spectral')
	for i in np.arange(0, x):
		for j in np.arange(0, y):
			new_cube[:,j,i] = rebin1d(new_wavels, wavels, datacube[:,j,i])
	
	
	#Update header
	head['CRVAL3'] = new_wavels[0]
	head['CDELT3'] = lamb_per_pix
	head['NAXIS3'] = len(new_wavels)
	head['SPECRES'] = new_res

	logging.info('# Spectral resolution and wavelength range - Done')
	
	#Return new_datacube, new_wavels, updated_header
	return new_cube, head, new_wavels, input_spec_res
	
	
def spatial_res(datacube, head, spax):
	'''Function that takes input datacube and rebins it to the
	chosen spatial resolution. 

	Inputs:
		datacube: object datacube
		head: header file for object datacube
		spax: spatial pixel (spaxel) scale 

	Outputs:
		new_cube: rebinned datacube
		head: updated header file
	'''

	logging.info('# Spatial resolution')
	logging.info('Chosen spaxel scale: ' + str(spax))

	z, y, x = datacube.shape

	try:
		spax_scale = config_data['spaxel_scale'][spax]
	except:
		raise HSIMError(spax + ' is not a valid spaxel scale. Valid options are: ' + ", ".join(sorted(config_data['spaxel_scale'].keys())))

	logging.info('Input sampling = %.2f x %.2f mas' % (head['CDELT1'], head['CDELT2']))
	new_sampling_x = spax_scale.psfscale*np.sign(head['CDELT1'])
	new_sampling_y = spax_scale.psfscale*np.sign(head['CDELT2'])
	
	x0 = head['CRVAL1'] + head['CDELT1']*(1 - head['CRPIX1'])
	y0 = head['CRVAL2'] + head['CDELT2']*(1 - head['CRPIX2'])
	
	logging.info('Internal sampling = %.2f x %.2f mas' % (new_sampling_x, new_sampling_y))

	if new_sampling_x != head['CDELT1'] or new_sampling_y != head['CDELT2']:
		
		# regrid image and conserve flux
		xmax = head['CDELT1']*(x-1)
		ymax = head['CDELT2']*(y-1)
			
		xgrid_in = np.linspace(0, xmax, x)
		ygrid_in = np.linspace(0, ymax, y)
			
		xgrid_out = np.arange(0, xmax, new_sampling_x)
		ygrid_out = np.arange(0, ymax, new_sampling_y)
		
		new_cube = np.zeros((z, len(ygrid_out), len(xgrid_out)), dtype=float)
		
		if new_sampling_x < head['CDELT1']:
			logging.info('Interpolating data cube - spatial')
			for k in np.arange(0, z):
				image = interp2d(xgrid_in, ygrid_in, datacube[k,:,:], kind='linear')
				new_cube[k,:,:] = image(xgrid_out, ygrid_out)
			
		else:
			logging.info('Rebinning data cube - spatial')
			for k in np.arange(0, z):
				new_cube[k,:,:] = frebin2d(datacube[k,:,:], (len(xgrid_out), len(ygrid_out)))
			
	else:
		
		new_cube = datacube
		
		
		
	#Update header
	head['CRPIX1'] = 1
	head['CRVAL1'] = x0
	head['CDELT1'] = new_sampling_x
	head['NAXIS1'] = new_cube.shape[2]
	
	head['CRPIX2'] = 1
	head['CRVAL2'] = y0
	head['CDELT2'] = new_sampling_y
	head['NAXIS2'] = new_cube.shape[1]
	

	logging.info('# Spatial resolution - Done')
	
	#Return new_datacube, updated_header
	return new_cube, head


def init_cube(datacube, grating, spax):
	''' Read input fits cube and resamples spatial and spectral 
		depending on selected grating and spaxel scale

	Inputs:
		datacube: Input high resolution datacube (RA, DEC, lambda)
		grating: Spectral grating
		spax: spatial pixel (spaxel) scale 

	Outputs:
		resampled cube

	'''
	
	#OPEN INPUT FITS file
	if os.path.isfile(datacube) == True and os.path.splitext(datacube)[1].lower() == '.fits':
		cube, head = fits.getdata(datacube, 0, header=True, memmap=True)
		#If not FITS file, try passing datacube directly
	else:
		try:
			cube = datacube.data
			head = datacube.header
		except:
			raise HSIMError('Please use FITS input file - ' + str(datacube))
	
	#Check that datacube has required headers to be processed in simulator
	#Required headers = ['CDELT1/2/3'], ['CRVAL3'], ['NAXIS1/2/3'], ['FUNITS'], ['CRPIX3'],
	#['CTYPE1/2/3'] = 'RA, DEC, WAVELENGTH, ['CUNIT1/2/3'] = MAS, MAS, microns/angstroms/etc,
	#['SPECRES'].
	fits_header_check(head)
	
	#CREATE WAVELENGTH ARRAY IN MICRONS
	lambs, head = wavelength_array(head)
	
	#RESCALE DATACUBE TO CHOSEN SPECTRAL RESOLUTION.
        cube, head, lambs, input_spec_res = spectral_res(cube, head, grating, lambs)
        
	#RESCALE DATACUBE TO CHOSEN SPATIAL RESOLUTION.
        cube, head = spatial_res(cube, head, spax)

	#Energy-to-Photons Conversion factor will depend on head['FUNITS'] value
	logging.info('Flux units = ' + head['FUNITS'])
	if head['FUNITS'] == 'J/s/m2/um/arcsec2':
		en2ph_conv_fac = (sp.h * sp.c)/(lambs*1.E-6) #J
	elif head['FUNITS'] == 'erg/s/cm2/A/arcsec2':
		en2ph_conv_fac = (sp.h * sp.c * 1.E7)/(lambs*1.E-6 * 1.E4 * 1.E4) #erg/1.E4(cm2->m2)/1.E4(A->um)
	elif head['FUNITS'] == 'J/s/m2/A/arcsec2':
		en2ph_conv_fac = (sp.h * sp.c)/(lambs*1.E-6 * 1.E4) #J/1.E4(A->um)
	elif head['FUNITS'] == 'erg/s/cm2/um/arcsec2':
		en2ph_conv_fac = (sp.h * sp.c * 1.E7)/(lambs*1.E-6 * 1.E4) #erg/1.E4(cm2->m2)
	else:
		raise HSIMError('Unknown flux units: Please change FUNITS header key to erg/s/cm2/A/arcsec2')
	
	en2ph_conv_fac.shape = (len(lambs),1,1)
	cube = np.divide(cube, en2ph_conv_fac)
	head['FUNITS'] = 'ph/s/m2/um/arcsec2'

	return cube, head, lambs, input_spec_res


