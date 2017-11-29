'''Code that contains the PSF functions used
in the HARMONI simulator.

Writen by Simon Zieleniewski

Started 04-06-13
Last edited 05-08-16
'''

import os
import numpy as n
import astropy.io.fits as p
import scipy.interpolate as s
from scipy.optimize import curve_fit
from modules.Gaussians import Gauss2D
from modules.frebin import *
from modules.spaxel_rebin import *
from modules.Math_functions import x1, x2, x6
from modules.misc_utils import path_setup
from modules.fits_utils import psf_fits_header_check, wavelength_array


psf_path = path_setup('../../Sim_data/PSFs/')


#AO options = 'SCAO', 'LTAO', 'Gaussian'
#If user provides PSF - this always supersedes AO choice
#PSF convolution ideally needs to happen at a spatial scale at least 1/10th of the output scale
#for the coarsest scales(10,20,40,60 mas) and at least 1/5th of the output scale of the finest
#scale (4, 5 mas).
#Depending on the input and output cube spatial scales, generate the PSF cube accordingly,
#perform the convolution with both input cube and PSF cube at same spatial scale then
#frebin up to the output spatial scale.

def psf_setup(cube, head, lambs, spax, user_PSF, AO, seeing, tel):
    '''Function to perform PSF setup

        Inputs:

             cube: datacube
             head: header
             spax: spaxel scale - tuple (x, y)
             user_PSF: path to 2D user PSF FITS file
             AO: AO mode: [LTAO, SCAO, Gaussian]
             seeing: seeing FWHM value [arcsec]
             zenith_ang: zenith angle [deg]
             tel: list of telescope diameter and obsc. ratio [D, eps]
             res_jitter: residual telescope windshake jitter [mas]
             Nyquist: False
             samp: PSF spatial sampling [mas]


        Outputs:

            cube: modified cube
            head: modified header
            psfspax: Initial PSF spaxel scale
            psfparams: dictionary containing several PSF parameters
            psfsize: array size in spaxels
            upsf: User uploaded PSF (2D or 3D)
            upsflams: Wavelength array for PSF cube if 3D cube


    '''

    #If user uploaded PSF
    if user_PSF != 'None':

        print 'User uploaded PSF'

        upsf, upsfh = p.getdata(user_PSF, header=True)

        psf_fits_header_check(upsfh)

        print 'Input PSF sampling = (%.2f, %.2f) mas' % (upsfh['CDELT1'], upsfh['CDELT2'])
        print 'Input datacube sampling = (%.2f, %.2f) mas' % (head['CDELT1'], head['CDELT2'])
        print 'Output spaxel scale = (%.2f, %.2f) mas' % (spax[0], spax[1])

        #Check spaxel scales of input and output cubes:
        mscale = 10.0
        if head['CDELT1'] <= spax[0]/mscale and head['CDELT2'] <= spax[1]/mscale:
            print 'Input spatial scale is at least 1/'+str(int(mscale))+'th of the output scale.'
            print 'Rebinning input cube to 1/'+str(int(mscale))+'th of output scale.'
            cube *= (head['CDELT1']*head['CDELT2']*1.E-6)
            if spax[0] <= spax[1]:
                print 'Output xspax < yspax'
                cube, head = spaxel_scale(cube, head, (spax[0]/mscale, spax[0]/mscale))
            elif spax[0] > spax[1]:
                print 'Output yspax < xspax'
                cube, head = spaxel_scale(cube, head, (spax[1]/mscale, spax[1]/mscale))
            cube /= (head['CDELT1']*head['CDELT2']*1.E-6)

        else:
            print 'WARNING: Input spatial scale is coarser than 1/'+str(int(mscale))+'th of the output scale.'
            print 'Interpolating input to a scale of 1/'+str(int(mscale))+'th the output.'
            print 'WARNING: Use of interpolation on input data!'
            cube *= (head['CDELT1']*head['CDELT2']*1.E-6)
            if spax[0] <= spax[1]:
                print 'Output xspax < yspax'
                cube, head = spaxel_scale(cube, head, (spax[0]/mscale, spax[0]/mscale))
            elif spax[0] > spax[1]:
                print 'Output yspax < xspax'
                cube, head = spaxel_scale(cube, head, (spax[1]/mscale, spax[1]/mscale))
            cube /= (head['CDELT1']*head['CDELT2']*1.E-6)

        print 'Scaling PSF to match input datacube'
        if upsfh['CDELT1'] > head['CDELT1'] and upsfh['CDELT2'] > head['CDELT2']:
            print 'PSF coarser than datacube - rebinning PSF to match'
            # cube *= (head['CDELT1']*head['CDELT2']*1.E-6)
            # cube, head = spaxel_scale(cube, head, (upsfh['CDELT1'],upsfh['CDELT2']))
            # cube /= (head['CDELT1']*head['CDELT2']*1.E-6)
            upsf, upsfh = spaxel_scale(upsf, upsfh, (head['CDELT1'],head['CDELT2']))
        if upsfh['CDELT1'] < head['CDELT1'] and upsfh['CDELT2'] < head['CDELT2']:
            print 'PSF finer than datacube - rebinning PSF to match'
            upsf, upsfh = spaxel_scale(upsf, upsfh, (head['CDELT1'],head['CDELT2']))

        #2D PSF image
        if 'NAXIS3' not in upsfh:
            upsflams = 'None'
        #3D PSF cube
        elif 'NAXIS3' in upsfh:
            upsflams, upsfh = wavelength_array(upsfh)
            if upsflams[0] <= lambs[0] and upsflams[-1] >= lambs[-1]:
                pass
            else:
                print 'Input PSF datacube DOES NOT cover wavelength range of input datacube'
                raise ValueError('Input PSF datacube DOES NOT cover wavelength range of input datacube')

        psfspax = upsfh['CDELT1']
        psfsize = upsfh['NAXIS1']
        psfparams = {} #placeholder


    else:
        upsf = 'None'; upsflams = 'None' #placeholders for when no user PSF given

        print 'Inbuilt PSFs'
        print 'Input spatial scales = %.1f mas, %.1f mas' % (head['CDELT1'], head['CDELT2'])
        print 'Chosen output spatial scales = %.1f mas, %.1f mas' % (spax[0], spax[1])
        psfspax = 1.0 #[mas/spaxel]
        print  'Chosen PSF initial sampling scale = %.1f mas' % psfspax

        #Check spaxel scales of input and output cubes:
        if head['CDELT1'] <= spax[0]/10. and head['CDELT2'] <= spax[1]/10.:
            print 'Input spatial scale is at least 1/10th of the output scale.'
            print 'Rebinning input cube to 1/10th of output scale.'
            cube *= (head['CDELT1']*head['CDELT2']*1.E-6)
            #cube, head = spaxel_scale(cube, head, (spax[0]/10., spax[1]/10.))
            if spax[0] <= spax[1]:
                print 'Output xspax < yspax'
                cube, head = spaxel_scale(cube, head, (spax[0]/10., spax[0]/10.))
            elif spax[0] > spax[1]:
                print 'Output yspax < xspax'
                cube, head = spaxel_scale(cube, head, (spax[1]/10., spax[1]/10.))
            cube /= (head['CDELT1']*head['CDELT2']*1.E-6)
            print 'Generating PSFcube at 1/10th output scale.'

            psfparams, psfsize = generate_psfcube(lambs, AO, seeing, tel, head, samp=psfspax)


        elif spax[0]/10. <= head['CDELT1'] <= spax[0]/5. and spax[1]/10. <= head['CDELT2'] <= spax[1]/5.:
            print 'Input spatial scale is between 1/10th and 1/5th of the output scale.'
            if spax[0] < 10. and spax[1] < 10.:
                print 'Output spaxel scale is finer than 10 mas. Convolution should be fine.'
                print 'Generating PSFcube at same spatial scale as input cube.'

                psfparams, psfsize = generate_psfcube(lambs, AO, seeing, tel, head, samp=psfspax)


            elif spax[0] >= 10. and spax[1] >= 10.:
                print 'Spaxel scale is equal to or coarser than 10 mas.'
                print 'WARNING: If coarser this runs the risk of convolving with the pixel scale twice!'
                print 'Generating PSFcube at same spatial scale as input cube.'

                psfparams, psfsize = generate_psfcube(lambs, AO, seeing, tel, head, samp=psfspax)

        else:
            print 'WARNING: Input spatial scale is coarser than 1/5th of the output scale.'
            print 'Interpolating input to a scale of 1/10th the output.'
            print 'WARNING: Use of interpolation on input data!'
            cube *= (head['CDELT1']*head['CDELT2']*1.E-6)
            #cube, head = spaxel_scale(cube, head, (spax[0]/10., spax[1]/10.))
            if spax[0] <= spax[1]:
                print 'Output xspax < yspax'
                cube, head = spaxel_scale(cube, head, (spax[0]/10., spax[0]/10.))
            elif spax[0] > spax[1]:
                print 'Output yspax < xspax'
                cube, head = spaxel_scale(cube, head, (spax[1]/10., spax[1]/10.))
            cube /= (head['CDELT1']*head['CDELT2']*1.E-6)

            psfparams, psfsize = generate_psfcube(lambs, AO, seeing, tel, head, samp=psfspax)

    print 'PSF setup done!'
    return cube, head, psfspax, psfparams, psfsize, upsf, upsflams



def generate_psfcube(lambs, AO, seeing, aperture, datahead, samp=1.):
    '''Function to generate PSF datacube parameters.

    Inputs:

        lambs: array of wavelengths for datacube
        AO: Type of AO: LTAO, SCAO, Gaussian. This option selects the
            parameters to generate chosen PSF type
        seeing: value seeing FWHM in arcsec.
        zenitha: zenith angle in degrees.
        aperture: List containing [diameter of telescope, obscuration ratio].
                  Diameter only used when Nyquist=True.
        res_jitter: Residual jitter corresponding to windshake/vibrations.
                    Entered in units of rms mas.
                    If Gaussian PSF is chosen, FWHM is given by
                    seeing FWHM + residual jitter.
        datahead: Datacube FITS header
        Nyquist: Boolean. Use nyquist sampling.
        samp: PSF sampling in milliarcsec.

    Outputs:

        psfcube: dictionary containing PSF parameters as a function of wavelength
        size: size of 2D PSF array in pixels
    '''

    print 'PSF parameter generation'


    #Get longest side of datacube array
    if datahead['NAXIS1'] >= datahead['NAXIS2']:
        datsize = datahead['NAXIS1']
        datsamp = datahead['CDELT1']
    else:
        datsize = datahead['NAXIS2']
        datsamp = datahead['CDELT2']


    if AO == 'LTAO':
        #Set array size according to: size = 2400 [mas] /spaxel_scale [mas/pixel] (10^-8 intensity achieved by r=1200 mas.)
        if datsize > int(2400/(datsamp*float(samp))):
            psfsize = datsize*datsamp/float(samp)
        else:
            #psfsize = int(2400/(datsamp*float(samp)))
            psfsize = 2400
        psfcube = LTAOpsfcube(lambs, seeing, aperture)

    elif AO == 'SCAO':
        #Set array size according to: size = 3000 [mas] /spaxel_scale [mas/pixel] (10^-7 intensity achieved by r=1200 mas.)
        if datsize > int(3000/(datsamp*float(samp))):
            psfsize = datsize*datsamp/float(samp)
        else:
            #psfsize = int(3000/(datsamp*float(samp)))
            psfsize = 3000
        psfcube = SCAOpsfcube(lambs, seeing, aperture)

    elif AO == 'Gaussian':
        #Set array size according to: size = 3200 [mas] /spaxel_scale [mas/pixel] (2x10^-8 intensity achieved by r=1200 mas.)
        if datsize > int(3200/(datsamp*float(samp))):
            psfsize = datsize*datsamp/float(samp)
        else:
            #psfsize = int(3200/(datsamp*float(samp)))
            psfsize = 3200
        psfcube = None

    else:
        return ValueError('Incorrect choice for AO type. Try again.')


    return psfcube, psfsize




def SCAOpsfcube(lambs, seeing, aperture):
    '''Function that creates a 3D spatial E-ELT SCAO PSF datacube by interpolating
    between parameters of several analytical functions. Parameters are stored in
    datafiles and are accessed by code.

    Inputs:

        lambs: Array of wavelengths.
        seeing: FWHM value of the seeing at 500nm (V-band). Value must be between 0.67" and 1.10".
        aperture: List containing [diameter of telescope, obscuration ratio]

    Output:

        psfcube: PSF cube of same spectral length as wavelength array
    '''

    print 'Generating SCAO PSF cube'

    ks = 1
    kl = 2
    box = [0.4,2.5,0.6,1.5]

    ###Load data from text file - wavelengths = vals[:,0]
    #vals = n.loadtxt('/Users/SimonZ/Data/Sim_data/PSFs/SCAOdata.txt', delimiter=',')
    #vals = n.loadtxt(psf_path+'SCAO/SCAOdata.txt', delimiter=',')
    vals = n.loadtxt(os.path.join(psf_path,'SCAO/SCAOdata.txt'), delimiter=',')

    #Seeing values
    see_vals = n.array([0.67, 0.85, 0.95, 1.10])

    params= []
    #Interpolate for all height parameters (oh = 1, 12, 23, 34; mh = 3, 14, 25, 36; lh = 6, 17, 28, 39; m2h = 9, 20, 31, 42)
    ohvals = n.array([vals[:,1], vals[:,12], vals[:,23], vals[:,34]])
    ohinterp = s.RectBivariateSpline(vals[:,0], see_vals, ohvals.transpose(),kx=kl, ky=ks, bbox=box)
    yoh = ohinterp(lambs, seeing)

    mhvals = n.array([vals[:,3], vals[:,14], vals[:,25], vals[:,36]])
    mhinterp = s.RectBivariateSpline(vals[:,0], see_vals, mhvals.transpose(),kx=kl, ky=ks, bbox=box)
    ymh = mhinterp(lambs, seeing)

    lhvals = n.array([vals[:,6], vals[:,17], vals[:,28], vals[:,39]])
    lhinterp = s.RectBivariateSpline(vals[:,0], see_vals, lhvals.transpose(),kx=kl, ky=ks, bbox=box)
    ylh = lhinterp(lambs, seeing)

    m2hvals = n.array([vals[:,9], vals[:,20], vals[:,31], vals[:,42]])
    m2hinterp = s.RectBivariateSpline(vals[:,0], see_vals, m2hvals.transpose(),kx=kl, ky=ks, bbox=box)
    ym2h = m2hinterp(lambs, seeing)

    params.append(yoh)
    params.append(ymh)
    params.append(ylh)
    params.append(ym2h)


    #Fit 6th order polynomial to Moffat (seeing) width and Lorentzian width parameters and then interpolate
    for i in n.array([4, 8]):

        p1, pn1 = curve_fit(x6, vals[:,0], vals[:,i])
        p2, pn2 = curve_fit(x6, vals[:,0], vals[:,i+11])
        p3, pn3 = curve_fit(x6, vals[:,0], vals[:,i+22])
        p4, pn4 = curve_fit(x6, vals[:,0], vals[:,i+33])
        yp1 = x6(vals[:,0], p1[0], p1[1], p1[2], p1[3], p1[4], p1[5], p1[6])
        yp2 = x6(vals[:,0], p2[0], p2[1], p2[2], p2[3], p2[4], p2[5], p2[6])
        yp3 = x6(vals[:,0], p3[0], p3[1], p3[2], p3[3], p3[4], p3[5], p3[6])
        yp4 = x6(vals[:,0], p4[0], p4[1], p4[2], p4[3], p4[4], p4[5], p4[6])

        yps = n.array([yp1, yp2, yp3, yp4])

        pinterp = s.RectBivariateSpline(vals[:,0], see_vals, yps.transpose(), kx=kl, ky=ks, bbox=box)
        ys = pinterp(lambs, seeing)

        params.append(ys)


    #Fit 1st order polynomial to Airy width, Moffat (seeing) shape, Lorentz position, Moffat(core) width
    for i in n.array([2, 5, 7, 10]):

        p1, pn1 = curve_fit(x1, vals[:,0], vals[:,i])
        p2, pn2 = curve_fit(x1, vals[:,0], vals[:,i+11])
        p3, pn3 = curve_fit(x1, vals[:,0], vals[:,i+22])
        p4, pn4 = curve_fit(x1, vals[:,0], vals[:,i+33])
        yp1 = x1(vals[:,0], p1[0], p1[1])
        yp2 = x1(vals[:,0], p2[0], p2[1])
        yp3 = x1(vals[:,0], p3[0], p3[1])
        yp4 = x1(vals[:,0], p4[0], p4[1])

        yps = n.array([yp1, yp2, yp3, yp4])
        pinterp = s.RectBivariateSpline(vals[:,0], see_vals, yps.transpose(), kx=kl, ky=ks, bbox=box)
        ys = pinterp(lambs, seeing)

        params.append(ys)

    #Fit 2nd order polynomial to Moffat (core) shape
    for i in n.array([11]):

        p1, pn1 = curve_fit(x2, vals[:,0], vals[:,i])
        p2, pn2 = curve_fit(x2, vals[:,0], vals[:,i+11])
        p3, pn3 = curve_fit(x2, vals[:,0], vals[:,i+22])
        p4, pn4 = curve_fit(x2, vals[:,0], vals[:,i+33])
        yp1 = x2(vals[:,0], p1[0], p1[1], p1[2])
        yp2 = x2(vals[:,0], p2[0], p2[1], p2[2])
        yp3 = x2(vals[:,0], p3[0], p3[1], p3[2])
        yp4 = x2(vals[:,0], p4[0], p4[1], p4[2])

        yps = n.array([yp1, yp2, yp3, yp4])
        pinterp = s.RectBivariateSpline(vals[:,0], see_vals, yps.transpose(), kx=kl, ky=ks, bbox=box)
        ys = pinterp(lambs, seeing)

        params.append(ys)

    #params = [oh, mh, lh, m2h, mw, lw, ow, mq, lp, m2w, m2q]
    pdict = {'oh':params[0], 'ow':params[6], 'mh':params[1],
             'mw':params[4], 'mq':params[7], 'lh':params[2],
             'lp':params[8], 'lw':params[5], 'm2h':params[3],
             'm2w':params[9], 'm2q':params[10]}
    return pdict




def LTAOpsfcube(lambs, seeing, aperture):
    '''Function that creates a 3D spatial E-ELT LTAO PSF datacube by interpolating
    between parameters of several analytical functions. Parameters are stored in
    datafiles and are accessed by code.

    Inputs:

        lambs: Array of wavelengths corresponding to [start_wave, end_wave, del_wave].
        seeing: FWHM value of the seeing. Value must be between 0.67" and 0.95" (21-12-13).
        aperture: List containing [diameter of telescope, obscuration ratio]

    Output:

        psfcube: PSF cube of same spectral length as wavelength array
    '''

    print 'Generating LTAO PSF cube'

    ks = 1
    kl = 2
    box = [0.4,2.5,0.6,1.5]

    ###Load data from text file - wavelengths = vals[:,0]
    #vals = n.loadtxt(psf_path+'LTAO/LTAOdata.txt', delimiter=',')
    vals = n.loadtxt(os.path.join(psf_path,'LTAO/LTAOdata.txt'), delimiter=',')

    #Seeing values
    see_vals = n.array([0.67, 0.95])

    params= []
    #Interpolate for all height parameters (oh = 1, 12, 23, 34; mh = 3, 14, 25, 36; lh = 6, 17, 28, 39; m2h = 9, 20, 31, 42)
    ohvals = n.array([vals[:,1], vals[:,11]])
    ohinterp = s.RectBivariateSpline(vals[:,0], see_vals, ohvals.transpose(),kx=kl, ky=ks, bbox=box)
    yoh = ohinterp(lambs, seeing)

    mhvals = n.array([vals[:,2], vals[:,12]])
    mhinterp = s.RectBivariateSpline(vals[:,0], see_vals, mhvals.transpose(),kx=kl, ky=ks, bbox=box)
    ymh = mhinterp(lambs, seeing)

    lhvals = n.array([vals[:,5], vals[:,15]])
    lhinterp = s.RectBivariateSpline(vals[:,0], see_vals, lhvals.transpose(),kx=kl, ky=ks, bbox=box)
    ylh = lhinterp(lambs, seeing)

    m2hvals = n.array([vals[:,8], vals[:,18]])
    m2hinterp = s.RectBivariateSpline(vals[:,0], see_vals, m2hvals.transpose(),kx=kl, ky=ks, bbox=box)
    ym2h = m2hinterp(lambs, seeing)

    params.append(yoh)
    params.append(ymh)
    params.append(ylh)
    params.append(ym2h)


    #Fit 6th order polynomial to Moffat (seeing) width and Lorentzian width parameters and then interpolate
    for i in n.array([3, 7]):

        p1, pn1 = curve_fit(x6, vals[:,0], vals[:,i])
        p2, pn2 = curve_fit(x6, vals[:,0], vals[:,i+10])
        yp1 = x6(vals[:,0], p1[0], p1[1], p1[2], p1[3], p1[4], p1[5], p1[6])
        yp2 = x6(vals[:,0], p2[0], p2[1], p2[2], p2[3], p2[4], p2[5], p2[6])

        yps = n.array([yp1, yp2])
        pinterp = s.RectBivariateSpline(vals[:,0], see_vals, yps.transpose(), kx=kl, ky=ks, bbox=box)
        ys = pinterp(lambs, seeing)

        params.append(ys)


    #Fit 1st order polynomial to Moffat (seeing) shape, Lorentz position
    for i in n.array([4, 6]):

        p1, pn1 = curve_fit(x1, vals[:,0], vals[:,i])
        p2, pn2 = curve_fit(x1, vals[:,0], vals[:,i+10])
        yp1 = x1(vals[:,0], p1[0], p1[1])
        yp2 = x1(vals[:,0], p2[0], p2[1])

        yps = n.array([yp1, yp2])
        pinterp = s.RectBivariateSpline(vals[:,0], see_vals, yps.transpose(), kx=kl, ky=ks, bbox=box)
        ys = pinterp(lambs, seeing)

        params.append(ys)

    #Interpolate for Moffat (core) width and shape as these are (effectively) step functions
    for i in n.array([9, 10]):
        yps = n.array([vals[:,i], vals[:,i+10]])
        pinterp = s.RectBivariateSpline(vals[:,0], see_vals, yps.transpose(), kx=kl, ky=ks, bbox=box)
        ys = pinterp(lambs, seeing)
        params.append(ys)

    #Airy width directly proportional to wavelength: width [mas] = (lambda[m]/(pi*D))*206265000.
    yows = lambs*1.E-6*206265000./(n.pi*aperture[0])
    params.append(yows)


    #params = [oh, mh, lh, m2h, mw, lw, mq, lp, m2w, m2q, ow]
    pdict = {'oh':params[0], 'ow':params[10], 'mh':params[1],
             'mw':params[4], 'mq':params[6], 'lh':params[2],
             'lp':params[7], 'lw':params[5], 'm2h':params[3],
             'm2w':params[8], 'm2q':params[9]}

    return pdict
