'''
Reads input FITS data cube and resamples according to the
grating and spaxel scale. The output cube is in ph/s/m2/um/arcsec2 units
'''
import os
import logging

import astropy.io.fits as fits
import astropy.units as u
import numpy as np
import scipy.constants as sp
from scipy.interpolate import interp2d

from src.config import *
from src.modules.rebin import *


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
		logging.warning("The input cube spectral resolution is lower than the HARMONI grating resolution. Assuming input resolution = HARMONI resolution")
		input_spec_res = new_res
	elif input_spec_res > 0.5*new_res:
		logging.warning("The input cube spectral resolution is lower than 2 times the HARMONI grating resolution")
	elif input_spec_res == 0.:
		logging.warning("The input cube spectral resolution is not defined")
	
	if input_spec_res < head['CDELT3']:
		logging.warning('Input resolution (%.1f AA) < Input sampling (%.1f AA). Assuming resultion = 2*sampling' % (input_spec_res*10000., head['CDELT3']*10000.))
		input_spec_res = 2.*head['CDELT3']
	
	logging.info('Input resolution = %.1f AA' % (input_spec_res*10000.))
	logging.info('Input sampling = %.1f AA' % (head['CDELT3']*10000.))
	logging.info('Output resolution = %.1f AA' % (new_res*10000.))
	logging.info('Internal sampling = %.1f AA' % (lamb_per_pix*10000.))

	#Interpolate datacube onto regular pixel grid
	if wavels[0] <= bandws.lmin and wavels[-1] >= bandws.lmax: # cube wavelength range larger than grating choice
		logging.warning('Cube wavelength range larger than grating: chopping down to grating range')
		new_wavels = np.arange(bandws.lmin, bandws.lmax, lamb_per_pix)
		
	elif wavels[0] > bandws.lmin and wavels[-1] < bandws.lmax: # cube wavelength range inside grating choice
		logging.info('Cube wavelength range within grating range')
		if len(wavels) > 1:
			new_wavels = np.arange(wavels[0], wavels[-1], lamb_per_pix)
		else:
			logging.warning('Input cube only has 1 slice. This slice will be duplicated in the internal cube')
			spec_size = 36
			new_wavels = np.arange(wavels[0] - spec_size*0.5*lamb_per_pix, wavels[0] + spec_size*0.5*lamb_per_pix, lamb_per_pix)
			wavels = new_wavels
	
	elif wavels[0] < bandws.lmin and wavels[-1] < bandws.lmax and wavels[-1] > bandws.lmin: # cube short wavelength longer than grating shortest wavelength
		logging.warning('Cube shortest wavelength shorter than grating')
		new_wavels = np.arange(bandws.lmin, wavels[-1], lamb_per_pix)
		
	elif wavels[0] > bandws.lmin and wavels[0] < bandws.lmax and wavels[-1] > bandws.lmax: # cube short wavelength longer than grating shortest wavelength
		logging.warning('Cube longest wavelength longer than grating')
		new_wavels = np.arange(wavels[0], bandws.lmax, lamb_per_pix)
		
	else:
		raise HSIMError('The wavelength range of the input cube ' + str(wavels[0]) + " - " + str(wavels[-1]) + " is not valid for the " + grating + " grating ("  + str(bandws.lmin) + " - " + str(bandws.lmax) + ")")

	
	# regrid spectrum and conserve flux
	if len(wavels) !=  len(new_wavels):
		logging.info('Interpolating data cube - spectral')
		new_cube = rebin_cube_1d(new_wavels, wavels, datacube)
	else:
		new_cube = datacube
	
	
	#Update header
	head['CRPIX3'] = 1
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

	min_internal_pix = 3.*spax_scale.yscale/spax_scale.psfscale
	
	xmax = head['CDELT1']*(x-1)
	ymax = head['CDELT2']*(y-1)
	if xmax < min_internal_pix or ymax < min_internal_pix:
		raise HSIMError('The input cube spatial dimension is too small. Minimum size is {s}x{s} mas'.format(s=int(min_internal_pix*spax_scale.psfscale)))
	

	if new_sampling_x != head['CDELT1'] or new_sampling_y != head['CDELT2']:
		
		# regrid image and conserve flux
		
		xgrid_in = np.linspace(0, xmax, x)
		ygrid_in = np.linspace(0, ymax, y)
		
		xgrid_out = np.arange(0, xmax, new_sampling_x)
		ygrid_out = np.arange(0, ymax, new_sampling_y)
		
		new_cube = np.zeros((z, len(ygrid_out), len(xgrid_out)), dtype=float)

		if new_sampling_x < head['CDELT1']:
			logging.warning('Interpolating data cube - spatial')
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
	
	if np.isnan(np.sum(cube)):
		raise HSIMError('NaN values are not allowed in the input cube')
	
	
	#Check that datacube has required headers to be processed in simulator
	required_headers = ['NAXIS1', 'NAXIS2', 'NAXIS3', 'CDELT1',
			'CDELT2', 'CDELT3', 'CRVAL3', 'BUNIT',
			'CRPIX3', 'CUNIT1', 'CUNIT2', 'CUNIT3',
			'CTYPE1', 'CTYPE2', 'CTYPE3', 'SPECRES']
	
	missing_headers = []
	
	for i in required_headers:
		if i not in head:
			logging.error('Missing header: ' + i)
			missing_headers.append(i)
			
	if len(missing_headers) != 0:
		raise HSIMError('Missing headers. Please correct datacube header.')
	
	# define CRVAL{1,2} and CRPIX{1,2} if not included in the header
	if "CRVAL1" not in head:
		head["CRVAL1"] = 0
	if "CRPIX1" not in head:
		head["CRPIX1"] = 1
		
	if "CRVAL2" not in head:
		head["CRVAL2"] = 0
	if "CRPIX2" not in head:
		head["CRPIX2"] = 1
	
	# Check axes types
	ctype1 = map(str.lower, ['ra', 'x', 'RA---SIN', 'RA---TAN'])
	ctype2 = map(str.lower, ['dec', 'y', 'DEC--SIN', 'DEC--TAN'])
	ctype3 = map(str.lower, ['wavelength'])
	
	if head['CTYPE1'].lower() not in ctype1:
		raise HSIMError("CTYPE1 must be set to any of the following: " + ", ".join(ctype1))
	
	if head['CTYPE2'].lower() not in ctype2:
		raise HSIMError("CTYPE2 must be set to any of the following: " + ", ".join(ctype2))
	
	if head['CTYPE3'].lower() not in ctype3:
		raise HSIMError("CTYPE3 must be set to any of the following: " + ", ".join(ctype3))
	
	# Check axes units
	try:
		head['CDELT1'] *= u.Unit(head["CUNIT1"]).to("mas")
		head['CUNIT1'] = 'mas'
	except ValueError as e:
		raise HSIMError("CUNIT1 error: " + str(e))
	
	try:
		head['CDELT2'] *= u.Unit(head["CUNIT2"]).to("mas")
		head['CUNIT2'] = 'mas'
	except ValueError as e:
		raise HSIMError("CUNIT2 error: " + str(e))
	
	try:
		factor = u.Unit(head["CUNIT3"]).to("micron")
		head['CDELT3'] *= factor
		head['SPECRES'] *= factor
		head['CRVAL3'] *= factor
		head['CUNIT3'] = 'micron'
	except ValueError as e:
		HSIMError("CUNIT3 error: " + str(e))
		
	
	# Create wavelength arrray 
	lambs = head['CRVAL3'] + head['CDELT3']*(np.linspace(1, head['NAXIS3'], head['NAXIS3']) - head['CRPIX3'])
	
	
	# Rescale datacube to chosen spectral resolution
	cube, head, lambs, input_spec_res = spectral_res(cube, head, grating, lambs)
	
	# Rescale datacube to chosen spatial resolution
	cube, head = spatial_res(cube, head, spax)


	z, y, x = cube.shape
	logging.info('Internal input cube size x={x} y={y} z={z}'.format(x=x, y=y, z=z))

	#Energy-to-Photons Conversion factor will depend on head['FUNITS'] value
	logging.info('Flux units = ' + head['BUNIT'])
	
	if head['BUNIT'] == "erg/s/cm2/A/arcsec2":
		head['BUNIT'] = "erg/s/cm2/AA/arcsec2"
		logging.warning("The input flux units should be erg/s/cm2/AA/arcsec2 instead of erg/s/cm2/A/arcsec2.")
	if head['BUNIT'] == "J/s/m2/A/arcsec2":
		head['BUNIT'] = "J/s/m2/AA/arcsec2"
		logging.warning("The input flux units should be J/s/m2/AA/arcsec2 instead of J/s/m2/A/arcsec2.")
		
	try:
		original_cube_units = u.Unit(head['BUNIT'])
		internal_cube_units = u.Unit("ph/s/m2/um/arcsec2")
		
		flux_factor = ((np.ones(len(lambs))*original_cube_units*u.arcsec**2).to(internal_cube_units*u.arcsec**2, equivalencies=u.spectral_density(lambs*u.micron))).base
		flux_factor.shape = (len(lambs),1,1)
		cube *= flux_factor
		
		head['BUNIT'] = str(internal_cube_units)
		
	except u.UnitConversionError as e:
		raise HSIMError("BUNIT error: " + str(e))
	
	logging.info('The flux range of the input cube is {:.2e} - {:.2e} ph/s/m2/um/arcsec2'.format(np.min(cube), np.max(cube)))

	spax_scale = config_data['spaxel_scale'][spax]
	
	area_spaxel = spax_scale.xscale*spax_scale.yscale/1000.**2 # arcsec2
	um_per_pixel = head['SPECRES'] # micron/pixel
	factor = config_data["telescope"]["area"]*area_spaxel*um_per_pixel
	
	logging.info('The flux range of the input cube is {:.2e} - {:.2e} ph/s/output pixel'.format(np.min(cube)*factor, np.max(cube)*factor))


	return cube, head, lambs, input_spec_res
