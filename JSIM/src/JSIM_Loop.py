'''Module that contains the loop over wavelength channels for the
simulation pipeline.

Written by Simon Zieleniewski

Created: 16-10-15

Last updated 02-12-15

'''

import numpy as n
from modules.frebin import *
from modules.Gaussians import Gauss2D
from modules.misc_utils import update_progress
from modules.psf_convolution import psf_convolve
import multiprocessing as mp
import webbpsf as wp

nspec = wp.NIRSpec()
#nspec.pupilopd = nspec.opd_list[0]
#nspec.pupilopd = None


def wavelength_loop(cube, head, wavels, out_cube, newsize, outspax):
    '''Function to take input datacube and process it iteratively through
    each wavelength channel as follows:
    - Generate PSF array for given channel
    - Add effect of ADR to channel
    - Convolve cube channel with PSF
    - Frebin up to chosen output spaxel scale

    Inputs:

        cube: Datacube
        head: Datacube header
        wavels: Wavelength array
        out_cube: Empty output cube
        newsize: tuple containing (x_newsize, y_newsize) array sizes
        outspax: output spaxels (x, y) (mas, mas)

    Output:

        cube: Processed cube
        head: Updated header

    '''
    nspec.pupilopd = None
    print 'OPD = ', nspec.pupilopd

    oversamp = 1000./float(outspax[0])

    for l in xrange(len(wavels)):
        #Print percentage bar
        update_progress(n.round(l/float(len(wavels)),2))

        #Create PSF channel
        psf = nspec.calcPSF(outfile=None, source=None, filter=None, nlambda=None,
                            monochromatic=wavels[l]*1.E-6, oversample=oversamp,
                            fov_arcsec=5, rebin=False)

        psf = psf[0].data

        #Convolve cube channel with PSF channel
        channel = psf_convolve(cube[l,:,:], psf)

        #Frebin datacube up to output spaxel scale
        newsize = (int(newsize[0]), int(newsize[1]))
        channel *= (head['CDELT1']*head['CDELT2']*1.E-6)
        channel = frebin(channel, (newsize[0],newsize[1]), total=True)
        channel /= (outspax[0]*outspax[1]*1.E-6)

        #Add channel to output cube
        out_cube[l,:,:] = channel

    return out_cube, head



def pp_wavelength_channel(chann, head, wavels, l, newsize, outspax):
    '''Function to take input datacube and process it iteratively through
    each wavelength channel as follows:
    - Generate PSF array for given channel
    - Add effect of ADR to channel
    - Convolve cube channel with PSF
    - Frebin up to chosen output spaxel scale

    Inputs:

        chann: cube channel
        head: Datacube header
        wave: wavelength [um]
        l: iteration no.
        out_cube: Empty output cube
        newsize: tuple containing (x_newsize, y_newsize) array sizes
        outspax: tuple containing (x, y) output spaxel scales

    Output:

        cube: Processed cube
        head: Updated header
        inspax: Input spaxel scale (mas, mas)
        outspax: Output spaxel scale (mas, mas)

    '''

    nspec.pupilopd = None
    print 'OPD = ', nspec.pupilopd

    oversamp = 1000./float(outspax[0])

    #Create PSF channel
    psf = nspec.calcPSF(outfile=None, source=None, filter=None, nlambda=None,
                        monochromatic=wavels[l]*1.E-6, oversample=oversamp,
                        fov_arcsec=5, rebin=False)
    psf = psf[0].data

    #Convolve cube channel with PSF channel
    chann = psf_convolve(chann, psf)

    #Frebin datacube up to output spaxel scale
    newsize = (int(newsize[0]), int(newsize[1]))
    chann *= (head['CDELT1']*head['CDELT2']*1.E-6)
    chann = frebin(chann, (newsize[0],newsize[1]), total=True)
    chann /= (outspax[0]*outspax[1]*1.E-6)

    return chann, l



def pp_wavelength_loop(cube, head, wavels, out_cube, newsize, outspax, usecpus=mp.cpu_count()-1):
    '''Function to take input datacube and process it iteratively through
    each wavelength channel as follows:
    - Generate PSF array for given channel
    - Add effect of ADR to channel
    - Convolve cube channel with PSF
    - Frebin up to chosen output spaxel scale

    Inputs:

        cube: Datacube
        head: Datacube header
        wavels: Wavelength array
        out_cube: Empty output cube
        newsize: tuple containing (x_newsize, y_newsize) array sizes
        outspax: tuple containing (x, y) output spaxel scales
        usecpus: no. of CPUs to use

    Output:

        cube: Processed cube
        head: Updated header
        inspax: Input spaxel scale (mas, mas)
        outspax: Output spaxel scale (mas, mas)

    '''
    import pprocess as pp

    print 'Using ', usecpus, ' CPUs'
    queue = pp.Queue(limit=usecpus)
    waveloop = queue.manage(pp.MakeParallel(pp_wavelength_channel))

    print "Calculating..."
    for lam in xrange(len(wavels)):
        #Print percentage bar
        #update_progress(n.round(lam/float(len(wavels)),2))
        waveloop(cube[lam,:,:], head, wavels[lam], lam, newsize, outspax)


    print "Finishing..."
    for chan, it in queue:
        update_progress(n.round(it/float(len(wavels)),2))
        out_cube[it,:,:] = chan


    return out_cube, head
