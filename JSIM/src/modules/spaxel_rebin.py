'''Spaxel rebinning function

Author: Simon Zieleniewski

Last updated: 19-05-15

'''

import numpy as n
from frebin import *


def spaxel_scale(datacube, head, spaxel):
    '''Function that takes datacube and rebins it to chosen
    spaxel scale. Reads spaxel scale of input datacube
    from header keywords CDELT1/2. 

    Inputs:
        datacube: science datacube
        head: header file
        spaxel: spaxel scale mas (x, y)

    Outputs:
        new_cube: datacube rebinned to chosen spaxel scale
        new_head: updated header file
    '''

    print 'Spaxel scale'

    z, y, x = datacube.shape
    cdelt1 = head["CDELT1"]
    cdelt2 = head["CDELT2"]

    #total field of view in mas
    x_field = cdelt1*x
    y_field = cdelt2*y

    x_newsize = n.round(x_field/float(spaxel[0]),0)
    y_newsize = n.round(y_field/float(spaxel[1]),0)

    newcube = n.zeros((z, y_newsize, x_newsize), dtype=n.float64)

    for i in xrange(z):
        newcube[i,:,:] = frebin(datacube[i,:,:], (x_newsize, y_newsize), total=True)

    head['CDELT1'] = (spaxel[0], "mas")
    head['CDELT2'] = (spaxel[1], "mas")
    head['NAXIS1'] = int(x_newsize)
    head['NAXIS2'] = int(y_newsize)

    print 'Spaxel scale - done!'
    
    return newcube, head


