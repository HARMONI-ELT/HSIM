'''FITS handling functions

Author: Simon Zieleniewski

Last updated: 11-02-16

'''

import numpy as n
import astropy.io.fits as p
import os


def generate_fits_cube(cube, head, wavels, filename, extra_headers, outd, varext=None):
    '''Function that takes a cube and generates a FITS file format
    version.

    Inputs:

        cube: datacube
        head: datacube header
        wavels: wavelength array
        filename: String of cube filename
        extra_headers: Dictionary containing extra header Keywords and values for FITS file.
        varext: Variance cube extension for Observed, Background and Reduced cubes

    Outputs:

        Saves FITS file in "Output_cubes" directory

    '''
    from misc_utils import path_setup

    #Update header with cube type
    #cube_type = filename.split('_')[0] + ' ' + filename.split('_')[1]
    cube_type = filename.split('_')[-2] + ' ' + filename.split('_')[-1]
    head['TYPE'] = cube_type

    #Update header with extras:
    head_keys = extra_headers.keys()
    for i in head_keys:
        head[i] = extra_headers[i]

    outFile = filename+'.fits'
    #os.chdir('./Output_cubes')
    cwd = os.getcwd()
    #os.chdir(path_setup('../../Output_cubes'))
    os.chdir(outd)

    for filename in os.listdir('./'):
        if (filename.lower() == outFile.lower()):
            os.rename(filename, filename.split('.fits')[0]+'.old'+'.fits')

    if head['CTYPE3'].lower() == 'wavelength':
        #p.writeto(outFile, cube, header=head)
        newfile = p.HDUList()
        newfile.append(p.PrimaryHDU(cube, head))
        if type(varext) == n.ndarray:
            vhead = head.copy()
            vhead['EXTTYPE'] = 'VARIANCE'
            newfile.append(p.ImageHDU(varext, vhead))
        newfile.writeto(outFile)

    elif head['CTYPE3'].lower() == 'nl wavelength':
        newfile = p.HDUList()
        newfile.append(p.PrimaryHDU(cube, head))
        whead = p.Header(['EXTTYPE', 'CUNIT1'])
        whead['EXTTYPE'] = 'WAVELENGTH'
        whead['CUNIT1'] = 'um'
        newfile.append(p.ImageHDU(wavels, whead))
        if type(varext) == n.ndarray:
            vhead = head.copy()
            vhead['EXTTYPE'] = 'VARIANCE'
            newfile.append(p.ImageHDU(varext, vhead))
        newfile.writeto(outFile)

    os.chdir(cwd)




def fits_header_check(header):
    '''Function that checks input FITS datacube header for
    required header keywords.
    Required headers = ['CDELT1/2/3'], ['CRVAL3'], ['NAXIS1/2/3'], ['FUNITS']/['BUNIT'], ['CRPIX3'] = 1,
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
    crpix3 = [1]

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
        print 'All required headers present!'

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
    if header['CRPIX3'] not in crpix3:
        raise HeaderError("CRPIX3 must be set to: ", crpix3)
    print 'All header values acceptable!'




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
    lambs = n.linspace(header['CRVAL3'], header['CRVAL3'] + header['CDELT3']*(header['NAXIS3']-1), header['NAXIS3'])

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
