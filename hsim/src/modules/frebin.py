'''Frebin function

Author: Simon Zieleniewski

Last updated: 19-05-15

'''

import numpy as n


def rebin(a, new_shape):
    """
    Resizes a 2d array by averaging or repeating elements, 
    new dimensions must be integral factors of original dimensions
 
    Parameters
    ----------
    a : array_like
        Input array.
    new_shape : tuple of int
        Shape of the output array (y, x)
 
    Returns
    -------
    rebinned_array : ndarray
        If the new shape is smaller of the input array, the data are averaged, 
        if the new shape is bigger array elements are repeated
 
    See Also
    --------
    resize : Return a new array with the specified shape.
 
    Examples
    --------
    >>> a = np.array([[0, 1], [2, 3]])
    >>> b = rebin(a, (4, 6)) #upsize
    >>> b
    array([[0, 0, 0, 1, 1, 1],
           [0, 0, 0, 1, 1, 1],
           [2, 2, 2, 3, 3, 3],
           [2, 2, 2, 3, 3, 3]])
 
    >>> c = rebin(b, (2, 3)) #downsize
    >>> c
    array([[ 0. ,  0.5,  1. ],
           [ 2. ,  2.5,  3. ]])
 
    """
    M, N = a.shape
    m, nn = new_shape
    if m<M:
        return a.reshape((m,M/m,nn,N/nn)).mean(3).mean(1)
    else:
        return n.repeat(n.repeat(a, m/M, axis=0), nn/N, axis=1)


    

def frebin(array, shape, total=True):
    '''Function that performs flux-conservative
    rebinning of an array.

    Inputs:
        array: numpy array to be rebinned
        shape: tuple (x,y) of new array size
	total: Boolean, when True flux is conserved

    Outputs:
	new_array: new rebinned array with dimensions: shape
    '''

    #Determine size of input image
    y, x = array.shape

    y1 = y-1
    x1 = x-1

    xbox = x/float(shape[0])
    ybox = y/float(shape[1])

    #Determine if integral contraction so we can use rebin
    if (shape[0] == int(shape[0])) and (shape[1] == int(shape[1])):
        if (x % shape[0] == 0) and (y % shape[1] == 0):
            return rebin(array, (shape[1], shape[0]))*xbox*ybox

    #Otherwise if not integral contraction
    #First bin in y dimension
    temp = n.zeros((int(shape[1]), x),dtype=n.float64)
    #Loop on output image lines
#    for i in xrange(0, int(n.round(shape[1],0)), 1):
    for i in xrange(0, int(shape[1]), 1):
        rstart = i*ybox
        istart = int(rstart)
        rstop = rstart + ybox
        istop = int(rstop)
        if istop > y1:
            istop = y1
        frac1 = rstart - istart
        frac2 = 1.0 - (rstop - istop)
        
    #Add pixel values from istart to istop an subtract
    #fracion pixel from istart to rstart and fraction
    #fraction pixel from rstop to istop.
        if istart == istop:
            temp[i,:] = (1.0 - frac1 - frac2)*array[istart,:]
        else:
            temp[i,:] = n.sum(array[istart:istop+1,:], axis=0)\
                        - frac1*array[istart,:]\
                        - frac2*array[istop,:]
            
    temp = n.transpose(temp)

    #Bin in x dimension
    result = n.zeros((int(shape[0]), int(shape[1])), dtype=n.float64)
    #Loop on output image samples
#    for i in xrange(0, int(n.round(shape[0],0)), 1):
    for i in xrange(0, int(shape[0]), 1):
        rstart = i*xbox
        istart = int(rstart)
        rstop = rstart + xbox
        istop = int(rstop)
        if istop > x1:
            istop = x1
        frac1 = rstart - istart
        frac2 = 1.0 - (rstop - istop)
    #Add pixel values from istart to istop an subtract
    #fracion pixel from istart to rstart and fraction
    #fraction pixel from rstop to istop.
        if istart == istop:
            result[i,:] = (1.-frac1-frac2)*temp[istart,:]
        else:
            result[i,:] = n.sum(temp[istart:istop+1,:], axis=0)\
                          - frac1*temp[istart,:]\
                          - frac2*temp[istop,:]

    if total:
        return n.transpose(result)
    elif not total:
        return n.transpose(result)/float(xbox*ybox)

