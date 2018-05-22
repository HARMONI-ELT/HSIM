'''
FITS handling functions
'''

import numpy as np

def fits_header_check(header):
	'''Function that checks input FITS datacube header for
	required header keywords.
	Required headers = ['CDELT1/2/3'], ['CDELT3'] > 0, ['CRVAL3'], ['NAXIS1/2/3'], ['FUNITS']/['BUNIT'], ['CRPIX3'],
	['CTYPE1/2/3'] = 'RA, DEC, WAVELENGTH, ['CUNIT1/2/3'] = MAS, MAS, microns/angstroms/etc,
	['SPECRES'].

		Input:

		header: FITS file header

	'''

	required_headers = ['NAXIS1', 'NAXIS2', 'NAXIS3', 'CDELT1',
				'CDELT2', 'CDELT3', 'CRVAL3', 'FUNITS',
				'CRPIX3', 'CUNIT1', 'CUNIT2', 'CUNIT3',
				'CTYPE1', 'CTYPE2', 'CTYPE3', 'SPECRES']
	missing_headers = []

	ctype1 = ['ra', 'x']
	ctype2 = ['dec', 'y']
	ctype3 = ['wavelength']
	cunit1 = ['mas', 'arcsec']
	cunit2 = cunit1
	cunit3 = ['micron', 'angstrom', 'meter', 'nanometers',
		'microns', 'angstroms', 'meters', 'nm']
	funits = ['J/s/m2/um/arcsec2', 'erg/s/cm2/A/arcsec2',
		'J/s/m2/A/arcsec2', 'erg/s/cm2/um/arcsec2']

	#Check for use of BUNIT key
	if 'BUNIT' in header:
		header['FUNITS'] = header['BUNIT']
		
	for i in required_headers:
		if i not in header:
			print 'Missing header: ', i
			missing_headers.append(i)
	
	if len(missing_headers) != 0:
		# print 'Missing headers: ', missing_headers
		raise HeaderError('Missing headers. Please correct datacube header.')
	else:
		if header['CUNIT1'].lower() == 'arcsec':
			header['CUNIT1'] = 'mas'
			header['CDELT1'] = header['CDELT1']*1000.
		if header['CUNIT2'].lower() == 'arcsec':
			header['CUNIT2'] = 'mas'
			header['CDELT2'] = header['CDELT2']*1000.
		print 'All required headers present'

	if header['CTYPE3'].lower() not in ctype3:
		raise HeaderError("CTYPE3 must be set to: ", ctype3)
	if header['CTYPE2'].lower() not in ctype2:
		raise HeaderError("CTYPE2 must be set to: ", ctype2)
	if header['CTYPE1'].lower() not in ctype1:
		raise HeaderError("CTYPE1 must be set to: ", ctype1)
	if header['CUNIT3'].lower() not in cunit3:
		raise HeaderError("CUNIT3 must be set to one of: ", cunit3)
	if header['CUNIT2'].lower() not in cunit2:
		raise HeaderError("CUNIT2 must be set to one of: ", cunit2)
	if header['CUNIT1'].lower() not in cunit1:
		raise HeaderError("CUNIT1 must be set to one of: ", cunit1)
	if header['FUNITS'] not in funits:
		raise HeaderError("FUNITS must be set to one of: ", funits)
	if header['CDELT3'] <= 0:
		raise HeaderError("CDELT3 must be positive")
	
	if "CRVAL1" not in header:
		header["CRVAL1"] = 0
	if "CRPIX1" not in header:
		header["CRPIX1"] = 1
		
	if "CRVAL2" not in header:
		header["CRVAL2"] = 0
	if "CRPIX2" not in header:
		header["CRPIX2"] = 1
	
	print 'All header values acceptable'


def wavelength_array(header):
	'''Function that creates wavelength array in microns
	using FITS file header.

		Inputs:

		header: FITS file header

		Outputs:

		lambs: wavelength array in microns
		head: updated FITS file header

	'''

	#Create wavelength array using FITS headers
	lambs = header['CRVAL3'] + header['CDELT3']*(np.linspace(1, header['NAXIS3'], header['NAXIS3']) - header['CRPIX3'])
	
	#Check wavelength value of headers ['CDELT3'], ['CRVAL3'], ['SPECRES'] using header ['CUNITS3']
	#Working units of simulator = MICRONS
	w_units = header['CUNIT3'].lower()

	if w_units == 'microns' or w_units == 'micron':
		pass
	elif w_units == 'angstroms' or w_units == 'angstrom':
		lambs *= 1.E-4
		header['SPECRES'] *= 1.E-4
	elif w_units == 'meters' or w_units == 'meter':
		lambs *= 1.E6
		header['SPECRES'] *= 1.E6
	elif w_units == 'nm' or w_units == 'nanometers':
		lambs *= 1.E-3
		header['SPECRES'] *= 1.E-3
	else:
		raise ValueError('Choose correct units please: microns, angstroms, metres or nm')
	#Update header with new wavelength units
	header['CRVAL3'] = lambs[0]
	header['CDELT3'] = (lambs[1]-lambs[0])
	header['CUNIT3'] = 'microns'

	return lambs, header



#Header Error function
class HeaderError(Exception):
	pass
