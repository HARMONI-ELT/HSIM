'''Generate cubes functions

Last updated: 03-09-15

'''

import numpy as n


#Object cube
def generate_object_cube(datacube, wavels, throughput_cube, DIT, NDIT, spaxel, delta_lambda, area, slow=False):
    '''Function that takes input datacube and other parameters and generates a cube of
    total object electrons per pixel.

    Inputs:
        datacube: object datacube
        wavels: array of wavelengths for datacube [um]
        throughput_cube: total throughput of system
        DIT: exposure time [s]
        NDIT: no. of exposures
        spaxel: tuple (x, y) spaxel sizes [arcsec].
        delta_lambda: spectral resolution element [um]
        area: telescope area [m2]

    Outputs:
        object_cube_new: cube representing total electrons per pixel for datacube
        object_noise: cube representing shot noise in electrons for each pixel
        object_cube: Noiseless object cube
    '''

    if slow==True:
        object_cube = datacube*throughput_cube*DIT*(spaxel[0]*spaxel[1])*delta_lambda*area

        object_cube_new = n.zeros(object_cube.shape, dtype=float)
        for x in xrange(NDIT):
            #Shot noise
            object_cube_new += n.random.poisson(object_cube)
        object_cube_new/=np.float(NDIT)
    
    else:
        object_cube = datacube*throughput_cube*DIT*NDIT*(spaxel[0]*spaxel[1])*delta_lambda*area
        object_cube_new = n.random.poisson(abs(object_cube))
        
    object_cube_new = object_cube_new.astype(n.float64)
    object_noise = n.abs(object_cube_new-object_cube)

    
    return object_cube_new, object_noise, object_cube


#Read cube
def generate_read_cube(cube_shape, wavels, DIT, NDIT, readsig_vis, readsig_nirs, slow=False):
    '''Function that creates read noise cube.
    Function uses read noise distribution from KMOS 300s
    dark frames.

    Inputs:

        cube_shape: (lambda,y,x)
        wavels: wavelength array [um]
        DIT: Exposure time [s]
        NDIT: no. of exposures
        readsig_vis: read noise visible [e/pix]
        readsig_nir: read noises for near-IR [high_exptime, low_exptime] [e/pix]

    Outputs:

        read_cube: Read noise cube
        
    '''

    #Read noise sigma from KMOS HAWAII-2RG detectors

    print 'Generating read noise cube'

    #Check what wavelengths cube is (Vis/near-IR cutoff at 0.8 um.)
    cutoff = 0.8

    #If DIT < 120 s use larger near-IR Read noise value
    if DIT <= 120.: 
        readsig_nir = readsig_nirs[1]
    else:
        readsig_nir = readsig_nirs[0]
    print 'DIT = ', DIT, 's, ', 'NIR READ = ', readsig_nir, 'e/pix'

    #Use normal distribution of mean=0., sigma=readsig
    #Also create "noiseless" read cube to use in SNR cubes
    read_cube = n.zeros(cube_shape, dtype=float)
    zero_cube = n.zeros_like(read_cube)
    nrc = n.zeros_like(read_cube)

    if slow==True:
        try:
            vis_cut = n.where(wavels < cutoff)[0][-1]
            for nnn in xrange(NDIT):
                read_cube[0:vis_cut+1,:,:] += n.random.normal(zero_cube[0:vis_cut+1,:,:], readsig_vis)
                nrc[0:vis_cut+1,:,:] += readsig_vis**2.0
                read_cube[vis_cut:,:,:] += n.random.normal(zero_cube[vis_cut:,:,:], readsig_nir)
                nrc[vis_cut:,:,:] += readsig_nir**2.0
            nrc=n.sqrt(nrc)
            
        except:
            #Check if wavelengths are shorter (visible) or longer (near-IR) than 0.8 um.
            if wavels[-1] < cutoff:
                for nnn in xrange(NDIT):
                    read_cube += n.random.normal(zero_cube, readsig_vis)
                    nrc += readsig_vis**2.0
            elif wavels[0] > cutoff:
                for nnn in xrange(NDIT):
                    read_cube += n.random.normal(zero_cube, readsig_nir)
                    nrc += readsig_nir**2.0
            nrc=n.sqrt(nrc)
                
    else:               
        try:
            vis_cut = n.where(wavels < cutoff)[0][-1]     
            read_cube[0:vis_cut+1,:,:] += n.random.normal(zero_cube[0:vis_cut+1,:,:], n.sqrt(NDIT)*readsig_vis)
            nrc[0:vis_cut+1,:,:] += n.sqrt(NDIT)*readsig_vis
            read_cube[vis_cut:,:,:] += n.random.normal(zero_cube[vis_cut:,:,:], n.sqrt(NDIT)*readsig_nir)
            nrc[vis_cut:,:,:] += n.sqrt(NDIT)*readsig_nir
        except:
            #Check if wavelengths are shorter (visible) or longer (near-IR) than 0.8 um.
            if wavels[-1] < cutoff:
                read_cube += n.random.normal(zero_cube, n.sqrt(NDIT)*readsig_vis)
                nrc += n.sqrt(NDIT)*readsig_vis
            elif wavels[0] > cutoff:
                read_cube += n.random.normal(zero_cube, n.sqrt(NDIT)*readsig_nir)
                nrc += n.sqrt(NDIT)*readsig_nir

    return read_cube, nrc



#Dark cube
def generate_dark_cube(cube_shape, wavels, dit, ndit, vis_dark, nir_dark, slow=False):
    '''Function that creates dark current cube.
    Currently only using placeholder value [e/s/pixel]
    from KMOS detectors.

    Inputs:

        cube_shape: (lambda, y, x)
        wavels: wavelength array [um]
        dit: Exposure time [s]
        ndit: No. of exposures

    Outputs:

        dark_cube_new: Cube of dark current [e-]
        dark_cube: Noiseless dark cube
    '''
    
    #Check what wavelengths cube is (Vis/near-IR cutoff at 0.8 um.)
    cutoff = 0.8

    dark_cube = n.zeros(cube_shape, dtype=float)
    try:
        vis_cut = n.where(wavels < cutoff)[0][-1]     
        dark_cube[0:vis_cut+1,:,:] += vis_dark
        dark_cube[vis_cut:,:,:] += nir_dark
    except:
        #Check if wavelengths are shorter (visible) or longer (near-IR) than 0.8 um.
        if wavels[-1] < cutoff:
            dark_cube += vis_dark
        elif wavels[0] > cutoff:
            dark_cube += nir_dark

    if slow==True:
        dark_cube *= dit

        dark_cube_new = n.zeros(dark_cube.shape, dtype=float)
        for x in xrange(ndit):
            #Shot noise
            dark_cube_new += n.random.poisson(dark_cube)
        dark_cube_new/=np.float(ndit)

    else:
        dark_cube *= dit*ndit
        dark_cube_new = n.random.poisson(abs(dark_cube))
        
    dark_cube_new = dark_cube_new.astype(n.float64)
    dark_noise = n.abs(dark_cube_new-dark_cube)

    return dark_cube_new, dark_cube



#Noise cube
def generate_noise_cube(obj_cube, bg_cube, dark_cube, read_cube):
    '''Function that takes input observed datacube and generates
    a cube corresponding to noise (sigma) for each pixel in the
    cube. Shot noise on object, background and dark current is
    computed using a standard normal distribution with sigma=1,
    mean=0. Read noise is an added constant dependent on
    wavelength and exposure time.

    Inputs:
        obj_cube: cube representing object variance
        bg_cube: cube representing background variance
        dark_cube: cube representing dark variance
        read_cube: cube representing read noise
        ndit: no. of detector integrations

    Outputs:
        noise_cube: cube respresnting sigma for each pixel in observed cube.
    '''
    
    #Uncorrelated noise adds in quadrature
    noise_cube = n.sqrt(obj_cube + bg_cube + dark_cube + read_cube**2)

    
    return noise_cube


