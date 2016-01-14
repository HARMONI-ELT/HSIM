'''Spatial PSF convolution function

Author: simon Zieleniewski

Last updated: 18-09-14

'''

import numpy as n



def psf_convolve(datacube, psfcube):
    '''Function to convolve each wavelength channel of datacube with
    corresponding spatial PSF.
    
    Inputs:
        datacube: science datacube
        psfcube: PSF datacube of same shape as science datacube
        
    Outputs:
        conv_cube: Convolved datacube of same shape as input science datacube
    '''

    yy, xx = psfcube.shape
    y, x = datacube.shape

    #If PSF images are equal or larger than cube images, pad cube images up to psf size
    if yy >= y and xx >= x:
#        print 'PSF cube equal or larger than datacube - padding datacube to PSF size - then FFT convolving!'
        extra_x = xx-x
        extra_y = yy-y
        xfac = 0; yfac = 0
        if extra_x % 2 != 0:
            extra_x += 1
            xfac = 1
        if extra_y % 2 != 0:
            extra_y += 1
            yfac = 1
        img_box = n.zeros((yy, xx), dtype=n.float64)

        img_box[(extra_y/2.)-yfac:yy-extra_y/2.,(extra_x/2.)-xfac:xx-extra_x/2.] = datacube
        conv_cube_slice = n.fft.ifft2(n.fft.fft2(img_box)*n.fft.fft2(psfcube))
        conv_cube_slice_final = n.fft.fftshift(conv_cube_slice.real)
        conv_channel = conv_cube_slice_final[(extra_y/2.)-yfac:yy-extra_y/2.,(extra_x/2.)-xfac:xx-extra_x/2.]
        
    #If PSF images are smaller than cube images, do something else!
    elif yy < y and xx < x:
#        print 'PSF cube smaller than datacube - padding PSF array with outermost values to PSF size - then FFT convolving.'
        extra_x = x-xx
        extra_y = y-yy
        xfac = 0; yfac = 0
        if extra_x % 2 != 0:
            extra_x += 1
            xfac = 1
        if extra_y % 2 != 0:
            extra_y += 1
            yfac = 1
        img_box = n.zeros((y, x), dtype=n.float64)

        img_box[(extra_y/2.)-yfac:y-extra_y/2.,(extra_x/2.)-xfac:x-extra_x/2.] = psfcube
        conv_cube_slice = n.fft.ifft2(n.fft.fft2(img_box)*n.fft.fft2(psfcube))
        conv_cube_slice_final = n.fft.fftshift(conv_cube_slice.real)
        conv_channel = conv_cube_slice_final[(extra_y/2.)-yfac:y-extra_y/2.,(extra_x/2.)-xfac:x-extra_x/2.]        

    else:
        print 'Datacube spatial shape = (%g, %g)' % (datacube.shape[1], datacube.shape[0])
        print 'PSF cube spatial shape = (%g, %g)' % (psfcube.shape[1], psfcube.shape[0])
        raise ValueError("Can't currently deal with odd shaped arrays!")
      
    return conv_channel

