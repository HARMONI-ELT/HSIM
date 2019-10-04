import numpy as n
import time
import sys

##sys.path.append('/Users/SimonZ/Documents/Oxford/D.Phil/Simulator/hsim1.0.0/src/modules/')
##import psf_convolution as psf
##import Gaussians as g
##sys.path.append('/Users/SimonZ/Documents/Oxford/D.Phil/Simulator/hsim1.0.0/src/cy_modules/src/cy_modules/')
##import cy_psf_convolution as cypsf


##x = g.Gauss2D(2000, 10)
##x.astype(n.double)
##y = n.zeros((2000,2000), dtype=n.double)
##y[y.shape[0]/2.,y.shape[1]/3.] = 1.
##y[y.shape[0]/2.,y.shape[1]*2./3.] = 1.
##
##
###Python
##print 'Python'
##start_time = time.time()
##z = psf.psf_convolve(x, y)
##print "Time elapsed: ", time.time() - start_time, "s"
##
##z = 0.
##
###Cython
##print 'Cython'
##start_time = time.time()
##cz = cypsf.psf_convolve(x, y)
##print "Time elapsed: ", time.time() - start_time, "s"
##
##cz = 0

##import frebin as f
##import cy_frebin as cf
##
##x = n.arange(0.,1000000., 1).reshape(1000,1000)
##
###Python
##print 'Python'
##start_time = time.time()
##z = f.frebin(x, n.array([325,325]))
##print "Time elapsed: ", time.time() - start_time, "s"
##
##z = 0.
##
###Cython
##print 'Cython'
##start_time = time.time()
##z = cf.frebin(x, n.array([325,32]))
##print "Time elapsed: ", time.time() - start_time, "s"


##sys.path.append('/Users/SimonZ/Documents/Oxford/D.Phil/Simulator/hsim1.0.0/src/')
##from TheSimulator_ADR import add_ADR
##import cy_add_ADR as cadr
##import astropy.io.fits as p
##import numpy as n
##
##cube, head = p.getdata('/Users/SimonZ/Documents/Oxford/D.Phil/Simulator/hsim1.0.0/HARMONI_Abs_line_gal_age01.0_Z0.02_logM10.0_z4.0_PointSource_v2.fits', header=True)
##a = 10.325
##
###Python
##print 'Python'
##start_time = time.time()
##z = add_ADR(cube[0,:,:], head, a)
##print "Time elapsed: ", time.time() - start_time, "s"
##
##z = 0.
##
###Cython
##print 'Cython'
##start_time = time.time()
##z = add_ADR(cube[0,:,:], head, a)
##print "Time elapsed: ", time.time() - start_time, "s"



##sys.path.append('/Users/SimonZ/Documents/Oxford/D.Phil/Simulator/hsim1.0.0/src/modules/')
##from psf_utils import create_Gausspsf_channel
##import numpy as n
##
##
###Python
##print 'Python'
##start_time = time.time()
##z = create_Gausspsf_channel(1.0, 0.7, n.array([37.,0.3]),
##                            5.0, 2000., False, 1.0, n.array([10.,10.]))
##print "Time elapsed: ", time.time() - start_time, "s"
##
##z = 0.
##
##sys.path.append('/Users/SimonZ/Documents/Oxford/D.Phil/Simulator/hsim1.0.0/src/cy_modules/src/cy_modules/')
##from cy_psf_utils import create_Gausspsf_channel
##
###Cython
##print 'Cython'
##start_time = time.time()
##z = create_Gausspsf_channel(1.0, 0.7, n.array([37.,0.3]),
##                            5.0, 2000., False, 1.0, n.array([10.,10.]))
##print "Time elapsed: ", time.time() - start_time, "s"


#Full wavelength loop speed test
if __name__=='__main__':
    import numpy as n
    import astropy.io.fits as p
    import sys
    sys.path.append('/Users/zieleniewski/Documents/Oxford/D.Phil/Simulator/hsim1.0.0/src/')
    import TheSimulator_Loop as sl
    import TheSimulator_PSFs as sp
    import TheSimulator_ADR as sa
    from TheSimulator_Config import *
    from src.modules.fits_utils import fits_header_check, wavelength_array

##    #Python
##    cube, head = p.getdata('/Users/zieleniewski/Documents/Oxford/D.Phil/Simulator/hsim1.0.0/HARMONI_Abs_line_gal_age01.0_Z0.02_logM10.0_z4.0_PointSource_v2.fits', header=True)
##    fits_header_check(head)
##    lambs, head = wavelength_array(head)
##    spax = n.array([100.,100.,])
##    out_cube = n.zeros((len(lambs),int(head['CDELT2']*head['NAXIS2']/spax[1]),
##                        int(head['CDELT1']*head['NAXIS1']/spax[0])), dtype=n.float64)
##    AO = 'LTAO'
##    seeing = 0.7
##    res_jitter = 0.0
##    user_PSF = 'None'
##    zenith_ang = 0.0
##    D = 37.
##    eps = 0.3
##    site_temp = 280.5
##    cube, head, psfspax, psfparams, psfsize, upsf = sp.psf_setup(cube, head, lambs, spax, user_PSF,
##                                                        AO, seeing, zenith_ang, n.array([D,eps]), res_jitter)
##    psfvars = [psfparams, psfspax, psfsize, n.array([D, eps]), res_jitter, seeing, upsf]
##    airmass = 1./n.cos(n.radians(zenith_ang))
##    refr = sa.calc_ref(lambs, site_temp)
##    adr = sa.calc_adr(lambs, refr, airmass)
##    print 'Python loop'
##    start_time = time.time()
##    out_cube, head = sl.wavelength_loop(cube, head, lambs, out_cube, AO, psfvars, adr,
##                                        (head['NAXIS1']*head['CDELT1']/float(spax[0]),head['NAXIS2']*head['CDELT2']/float(spax[1])),
##                                        spax, adr_switch=config_data['ADR'])
##    print "Time elapsed: ", time.time() - start_time, "s"
##    cube = 0.
##    out_cube = 0.


    #Parallel Python
    cube, head = p.getdata('/Users/zieleniewski/Documents/Oxford/D.Phil/Simulator/hsim1.0.0/HARMONI_Abs_line_gal_age01.0_Z0.02_logM10.0_z4.0_PointSource_v2.fits', header=True)
    fits_header_check(head)
    lambs, head = wavelength_array(head)
    spax = n.array([100.,100.,])
    out_cube = n.zeros((len(lambs),int(head['CDELT2']*head['NAXIS2']/spax[1]),
                        int(head['CDELT1']*head['NAXIS1']/spax[0])), dtype=n.float64)
    AO = 'LTAO'
    seeing = 0.7
    res_jitter = 0.0
    user_PSF = 'None'
    zenith_ang = 0.0
    D = 37.
    eps = 0.3
    site_temp = 280.5
    cube, head, psfspax, psfparams, psfsize, upsf = sp.psf_setup(cube, head, lambs, spax, user_PSF,
                                                        AO, seeing, zenith_ang, n.array([D,eps]), res_jitter)
    psfvars = [psfparams, psfspax, psfsize, n.array([D, eps]), res_jitter, seeing, upsf]
    airmass = 1./n.cos(n.radians(zenith_ang))
    refr = sa.calc_ref(lambs, site_temp)
    adr = sa.calc_adr(lambs, refr, airmass)
    print 'Parallel Python loop'
    start_time = time.time()
    out_cube, head = sl.pp_wavelength_loop(cube, head, lambs, out_cube, AO, psfvars, adr,
                                        (head['NAXIS1']*head['CDELT1']/float(spax[0]),head['NAXIS2']*head['CDELT2']/float(spax[1])),
                                        spax, adr_switch=config_data['ADR'])
    print "Time elapsed: ", time.time() - start_time, "s"
    cube = 0.
    out_cube = 0.
    

##    #Cython
##    sys.path.append('/Users/zieleniewski/Documents/Oxford/D.Phil/Simulator/hsim1.0.0/src/cy_modules/src/cy_modules')
##    from cy_wavelength_loop import wavelength_loop
##    cube, head = p.getdata('/Users/zieleniewski/Documents/Oxford/D.Phil/Simulator/hsim1.0.0/HARMONI_Abs_line_gal_age01.0_Z0.02_logM10.0_z4.0_PointSource_v2.fits', header=True)
##    fits_header_check(head)
##    lambs, head = wavelength_array(head)
##    spax = n.array([100.,100.,])
##    out_cube = n.zeros((len(lambs),int(head['CDELT2']*head['NAXIS2']/spax[1]),
##                        int(head['CDELT1']*head['NAXIS1']/spax[0])), dtype=n.float64)
##    AO = 'LTAO'
##    seeing = 0.7
##    res_jitter = 0.0
##    user_PSF = 'None'
##    zenith_ang = 0.0
##    D = 37.
##    eps = 0.3
##    site_temp = 280.5
##    cube, head, psfspax, psfparams, psfsize, upsf = sp.psf_setup(cube, head, lambs, spax, user_PSF,
##                                                        AO, seeing, zenith_ang, n.array([D,eps]), res_jitter)
##    psfvars = [psfparams, psfspax, psfsize, n.array([D, eps]), res_jitter, seeing, upsf]
##    airmass = 1./n.cos(n.radians(zenith_ang))
##    refr = sa.calc_ref(lambs, site_temp)
##    adr = sa.calc_adr(lambs, refr, airmass)
##    newsize = n.array([int(head['NAXIS1']*head['CDELT1']/float(spax[0])),int(head['NAXIS2']*head['CDELT2']/float(spax[1]))])
##    print 'Cython loop'
##    start_time = time.time()
##    out_cube, head = wavelength_loop(cube, head, lambs, out_cube, AO, psfvars, adr,
##                                     newsize, spax, adr_switch=config_data['ADR'])
##    print "Time elapsed: ", time.time() - start_time, "s"







