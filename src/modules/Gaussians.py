'''Gaussian functions

Author: Simon Zieleniewski

Last updated: 17-09-14

'''

import numpy as n

#1D
def Gauss(sigma, delta_x, limit=5):
    '''Function that creates an array of a normalised
    Gaussian distribution for use in convolutions.

    Inputs:
    
        sigma: Gaussian dispersion (related to FWHM by
               FWHM = 2*sqrt(2ln(2))*sigma
        delta_x: resolution [x/pixel]
        limit: Gaussian extent = [-limit*sigma, +limit*sigma]. Default=5

    Outputs:
    
        column array of (x, y)
    '''
    if sigma != 0.:
        #Create Gaussian
        num_xs = n.round(sigma*10/float(delta_x),decimals=0)
        if n.mod(num_xs,2) == 0.:
            num_xs += 1.
        xs = n.linspace(-5*sigma, 5*sigma, num=num_xs)
        ys = (1/(n.sqrt(2.*n.pi)*sigma))*n.exp(-0.5*(xs/float(sigma))**2)

        #normalise
        ys /= n.sum(ys)

        return n.column_stack((xs, ys))
    elif sigma == 0.:
        return n.array([1.,1.])


#2D
def Gauss2D(size, fwhm=3, center=None):
    ''' Creates a square normalised Gaussian kernel.

    Inputs:

        size: length of a side of the square grid
        fwhm: full-width-half-maximum value
        centre: (x,y) position of central Gaussian peak

    Outputs:

        arr: 2D array of normalised Gaussian
    '''

    x = n.arange(0, size, 1, n.float64)
    y = x[:,n.newaxis]

    if center is None:
        x0 = y0 = size // 2
    else:
        x0 = center[0]
        y0 = center[1]

    arr =  n.exp(-4*n.log(2) * ((x-x0)**2 + (y-y0)**2) / fwhm**2)
    arr /= n.sum(arr)
    
    return arr
