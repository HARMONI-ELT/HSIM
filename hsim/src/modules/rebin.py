'''
Rebin 1d and 2d arrays
'''


import numpy as np

from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
from scipy.integrate import quad

def rebin1d(xout, xin, yin):
	in0 = int(np.interp(xout[0], xin, range(len(xin))))
	
	dx_in = xin[in0+1] - xin[in0]
	dx_out = xout[1] - xout[0]

	if dx_out < dx_in:
		# interpolate if output is finer
		return np.interp(xout, xin, yin)
	else:
		# rebin if output is coarser
		temp = np.zeros((len(xout)), dtype=np.float64)
		#Loop on output values
		box = float(dx_out)/float(dx_in)
		
		in_i = np.interp(xout - dx_out*0.5, xin, range(len(xin)))
		
		for i in range(len(xout)):
			rstart = in_i[i]
			istart = int(rstart)
			if i < len(xout) - 1:
				rstop = in_i[i+1]
			else:
				# for the last one assume the previous box size
				rstop = in_i[i] + (in_i[i] - in_i[i-1])
				
			istop = int(rstop)
			if istop > len(xin) - 1:
				istop = len(xin) - 1
				
			frac1 = rstart - istart
			frac2 = 1.0 - (rstop - istop)
			
			#print istart, istop, rstart, rstop
			#Add pixel values from istart to istop an subtract
			#fracion pixel from istart to rstart and fraction
			#fraction pixel from rstop to istop.
			if istart == istop:
				temp[i] = (1.0 - frac1 - frac2)*yin[istart]
			else:
				temp[i] = np.sum(yin[istart:istop+1]) - frac1*yin[istart] - frac2*yin[istop]
		
		return np.transpose(temp)/box



def frebin2d(array, shape):
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


	#Otherwise if not integral contraction
	#First bin in y dimension
	temp = np.zeros((int(shape[1]), x), dtype=np.float64)
	#Loop on output image lines
	#    for i in xrange(0, int(np.round(shape[1],0)), 1):
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
			temp[i,:] = np.sum(array[istart:istop+1,:], axis=0)\
				- frac1*array[istart,:]\
				- frac2*array[istop,:]
		
	temp = np.transpose(temp)

	#Bin in x dimension
	result = np.zeros((int(shape[0]), int(shape[1])), dtype=np.float64)
	#Loop on output image samples
	#    for i in xrange(0, int(np.round(shape[0],0)), 1):
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
			result[i,:] = np.sum(temp[istart:istop+1,:], axis=0)\
				- frac1*temp[istart,:]\
				- frac2*temp[istop,:]

	return np.transpose(result)/float(xbox*ybox)


