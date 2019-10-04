"""
NGHXRG - Teledyne HxRG Noise Generator

Modification History:

8-15 April 2015, B.J. Rauscher, NASA/GSFC
- Implement Pierre Ferruit's (ESA/ESTEC)recommendation to use
  numpy.fft.rfft and numpy.fft.irfft for faster execution. This saves
  about 30% in execution time.
- Clean up setting default values per several suggestions from
  Chaz Shapiro (NASA/JPL).
- Change from setting the shell variable H2RG_PCA0 to point to the PCA-zero
  file to having a shell variable point to the NG home directory. This is per
  a suggestion from Chaz Shapiro (NASA/JPL) that would allow seamlessly adding  
  more PCA-zero templates.
- Implement a request form Chaz Shapiro (NASA/JPL) for status
  reporting. This was done by adding the "verbose" arguement.
- Implement a request from Pierre Ferruit (ESA/ESTEC) to generate
  3-dimensional data cubes.
- Implement a request from Pierre Ferruit to treat ACN as different 1/f
  noise in even/odd columns. Previously ACN was treated purely as a feature
  in Fourier space.
- Version 2(Beta)

16 April 2015, B.J. Rauscher
- Fixed a bug in the pinkening filter definitions. Abs() was used where
  sqrt() was intended. The bug caused power spectra to have the wrong shape at
  low frequency.
- Version 2.1(Beta)

17 April 2015, B.J. Rauscher
- Implement a request from Chaz Shapiro for HXRGNoise() to exit gracefully if
  the pca0_file is not found.
- Version 2.2 (Beta)

8 July 2015, B.J. Rauscher
- Address PASP referee comments
    * Fast scan direction is now reversible. To reverse the slow scan
      direction use the numpy flipud() function.
    * Modifications to support subarrays. Specifically,
        > Setting reference_pixel_border_width=0 (integer zero);
            + (1) eliminates the reference pixel border and
            + (2) turns off adding in a bias pattern
                  when simulating data cubes. (Turned on in v2.5)
- Version 2.4

12 Oct 2015, J.M. Leisenring, UA/Steward
- Make compatible with Python 2.x
    * from __future__ import division
- Included options for subarray modes (FULL, WINDOW, and STRIPE)
    * Keywords x0 and y0 define a subarray position (lower left corner)
    * Selects correct pca0 region for subarray underlay
    * Adds reference pixels if they exist within subarray window
- Tie negative values to 0 and anything >=2^16 to 2^16-1
- Version 2.5

20 Oct 2015, J.M. Leisenring, UA/Steward
- Padded nstep to the next power of 2 in order to improve FFT runtime
    * nstep2 = int(2**np.ceil(np.log2(nstep)))
    * Speeds up FFT calculations by ~5x
- Don't generate noise elements if their magnitudes are equal to 0.
- mknoise() now returns copy of final HDU result for easy retrieval
- Version 2.6

"""
# Necessary for Python 2.6 and later (JML)
# Should still work under Python 3.x (JML)
from __future__ import division, print_function

import os
import warnings
from astropy.io import fits
import numpy as np
from scipy.ndimage.interpolation import zoom
from scipy.interpolate import UnivariateSpline
from astropy.stats.funcs import median_absolute_deviation as mad
import datetime
# import matplotlib.pyplot as plt # Handy for debugging

from src.modules.misc_utils import path_setup
from src.config import *

warnings.filterwarnings('ignore')

# Have not verified this in Python 3.x (JML)
import logging
_log = logging.getLogger('nghxrg')

class HXRGNoise:
	"""
	HXRGNoise is a class for making realistic Teledyne HxRG system
	noise. The noise model includes correlated, uncorrelated,
	stationary, and non-stationary components. The default parameters
	make noise that resembles Channel 1 of JWST NIRSpec. NIRSpec uses
	H2RG detectors. They are read out using four video outputs at
	1.e+5 pix/s/output.
	"""

	# These class variables are common to all HxRG detectors
	nghxrg_version = 2.7 # Sofware version

	def __init__(self, naxis1=None, naxis2=None, naxis3=None, n_out=None,
				 dt=None, nroh=None, nfoh=None, pca0_file=None, verbose=False,
				 reverse_scan_direction=False, reference_pixel_border_width=None,
				 wind_mode='FULL', x0=0, y0=0, det_size=None, rn_file=None):
		"""
		Simulate Teledyne HxRG+SIDECAR ASIC system noise.

		Parameters:
			naxis1      - X-dimension of the FITS cube
			naxis2      - Y-dimension of the FITS cube
			naxis3      - Z-dimension of the FITS cube
						  (number of up-the-ramp samples)
			n_out       - Number of detector outputs
			nfoh        - New frame overhead in rows. This allows for a short
						  wait at the end of a frame before starting the next
						  one.
			nroh        - New row overhead in pixels. This allows for a short
						  wait at the end of a row before starting the next one.
			dt          - Pixel dwell time in seconds
			pca0_file   - Name of a FITS file that contains PCA-zero
			rn_file     - Name of a FITS file that contains read noisepca
			verbose     - Enable this to provide status reporting
			wind_mode   - 'FULL', 'STRIPE', or 'WINDOW' (JML)
			x0/y0       - Pixel positions of subarray mode (JML)
			det_size    - Pixel dimension of full detector (square), used only
						  for WINDOW mode (JML)
			reference_pixel_border_width - Width of reference pixel border
										   around image area
			reverse_scan_direction - Enable this to reverse the fast scanner
									 readout directions. This
									 capability was added to support
									 Teledyne's programmable fast scan
									 readout directions. The default
									 setting =False corresponds to
									 what HxRG detectors default to
									 upon power up.
		"""

		# ======================================================================
		#
		# DEFAULT CLOCKING PARAMETERS
		#
		# The following parameters define the default HxRG clocking pattern. The
		# parameters that define the default noise model are defined in the
		# mknoise() method.
		#
		# ======================================================================

		# Subarray Mode? (JML)
		if wind_mode is None:
			wind_mode = 'FULL'
		if det_size is None:
			det_size = 4096
		wind_mode = wind_mode.upper()
		modes = ['FULL', 'STRIPE', 'WINDOW']
		if wind_mode not in modes:
			_log.warn('%s not a valid window readout mode! Returning...' % inst_params['wind_mode'])
			os.sys.exit()
		if wind_mode == 'WINDOW':
			n_out = 1
		if wind_mode == 'FULL':
			x0 = 0; y0 = 0
		if wind_mode == 'STRIPE':
			x0 = 0

		# Default clocking pattern is JWST NIRSpec
		self.naxis1    = 4096  if naxis1   is None else naxis1
		self.naxis2    = 4096  if naxis2   is None else naxis2
		self.naxis3    = 1     if naxis3   is None else naxis3
		self.n_out     = 64     if n_out    is None else n_out
		self.dt        = 1.e-5 if dt       is None else dt
		self.nroh      = 12    if nroh     is None else nroh
		self.nfoh      = 1     if nfoh     is None else nfoh
		self.reference_pixel_border_width = 4 if reference_pixel_border_width is None \
											  else reference_pixel_border_width
											  
		# Check that det_size is greater than self.naxis1 and self.naxis2 in WINDOW mode (JML)
		if wind_mode == 'WINDOW':
			if (self.naxis1 > det_size):
				_log.warn('NAXIS1 %s greater than det_size %s! Returning...' % (self.naxis1,det_size))
				os.sys.exit()
			if (self.naxis2 > det_size):
				_log.warn('NAXIS2 %s greater than det_size %s! Returning...' % (self.naxis1,det_size))
				os.sys.exit()
				
                NGHXRG_HOME = path_setup('../../' + config_data["data_dir"] + 'detectors/')

		# Initialize PCA-zero file and make sure that it exists and is a file
		#self.pca0_file = os.getenv('NGHXRG_HOME')+'/nirspec_pca0.fits' if \
		#	pca0_file is None else pca0_file
		self.pca0_file = NGHXRG_HOME+'/nirspec_pca0.fits' if \
                        pca0_file is None else pca0_file
		if os.path.isfile(self.pca0_file) is False:
			print('There was an error finding pca0_file! Check to be')
			print('sure that the NGHXRG_HOME shell environment')
			print('variable is set correctly and that the')
			print('$NGHXRG_HOME/ directory contains the desired PCA0')
			print('file. The default is nirspec_pca0.fits.')
			os.sys.exit()

		# Initialize read noise file and make sure that it exists and is a file
		#self.rn_file = os.getenv('NGHXRG_HOME')+'/kmos_rn.fits' if \
                #        rn_file is None else rn_file
                self.rn_file = NGHXRG_HOME+'/kmos_rn.fits' if \
                        rn_file is None else rn_file
		if os.path.isfile(self.rn_file) is False:
                        print('There was an error finding rn_file! Check to be')
                        print('sure that the NGHXRG_HOME shell environment')
                        print('variable is set correctly and that the')
                        print('$NGHXRG_HOME/ directory contains the desired read')
                        print('noise file. The default is kmos_rn.fits.')
                        os.sys.exit()


		# ======================================================================

		# Configure Subarray (JML)
		self.wind_mode = wind_mode
		self.det_size  = det_size
		self.x0 = x0
		self.y0 = y0

		# Configure status reporting
		self.verbose = verbose

		# Configure readout direction
		self.reverse_scan_direction = reverse_scan_direction
	
		# Compute the number of pixels in the fast-scan direction per
		# output
		self.xsize = self.naxis1 // self.n_out
			
		# Compute the number of time steps per integration, per
		# output
		self.nstep = (self.xsize+self.nroh) * (self.naxis2+self.nfoh) * self.naxis3
		# Pad nsteps to a power of 2, which is much faster (JML)
		self.nstep2 = int(2**np.ceil(np.log2(self.nstep)))

		# For adding in ACN, it is handy to have masks of the even
		# and odd pixels on one output neglecting any gaps
		self.m_even = np.zeros((self.naxis3,self.naxis2,self.xsize))
		self.m_odd = np.zeros_like(self.m_even)
		for x in np.arange(0,self.xsize,2):
			self.m_even[:,:self.naxis2,x] = 1
			self.m_odd[:,:self.naxis2,x+1] = 1
		self.m_even = np.reshape(self.m_even, np.size(self.m_even))
		self.m_odd = np.reshape(self.m_odd, np.size(self.m_odd))

		# Also for adding in ACN, we need a mask that point to just
		# the real pixels in ordered vectors of just the even or odd
		# pixels
		self.m_short = np.zeros((self.naxis3, self.naxis2+self.nfoh, \
									  (self.xsize+self.nroh)//2))
		self.m_short[:,:self.naxis2,:self.xsize//2] = 1
		self.m_short = np.reshape(self.m_short, np.size(self.m_short))

		# Define frequency arrays
		self.f1 = np.fft.rfftfreq(self.nstep2) # Frequencies for nstep elements
		self.f2 = np.fft.rfftfreq(2*self.nstep2) # ... for 2*nstep elements

		# Define pinkening filters. F1 and p_filter1 are used to
		# generate ACN. F2 and p_filter2 are used to generate 1/f noise.
		self.alpha = -1 # Hard code for 1/f noise until proven otherwise
		self.p_filter1 = np.sqrt(self.f1**self.alpha)
		self.p_filter2 = np.sqrt(self.f2**self.alpha)
		self.p_filter1[0] = 0.
		self.p_filter2[0] = 0.


		# Initialize pca0. This includes scaling to the correct size,
		# zero offsetting, and renormalization. We use robust statistics
		# because pca0 is real data
		hdu = fits.open(self.pca0_file)
		nx_pca0 = hdu[0].header['naxis1']
		ny_pca0 = hdu[0].header['naxis2']

		# Do this slightly differently, taking into account the 
		# different types of readout modes (JML)
		#if (nx_pca0 != self.naxis1 or naxis2 != self.naxis2):
		#    zoom_factor = self.naxis1 / nx_pca0
		#    self.pca0 = zoom(hdu[0].data, zoom_factor, order=1, mode='wrap')
		#else:
		#    self.pca0 = hdu[0].data
		#self.pca0 -= np.median(self.pca0) # Zero offset
		#self.pca0 /= (1.4826*mad(self.pca0)) # Renormalize

		data = hdu[0].data        
		# Make sure the real PCA image is correctly scaled to size of fake data (JML)
		# Depends if we're FULL, STRIPE, or WINDOW
		if wind_mode == 'FULL':
			scale1 = self.naxis1 / nx_pca0
			scale2 = self.naxis2 / ny_pca0
			zoom_factor = np.max([scale1, scale2])
		if wind_mode == 'STRIPE':
			zoom_factor = self.naxis1 / nx_pca0
		if wind_mode == 'WINDOW':
			# Scale based on det_size
			scale1 = self.det_size / nx_pca0
			scale2 = self.det_size / ny_pca0
			zoom_factor = np.max([scale1, scale2])
		
		# Resize PCA0 data
		if zoom_factor != 1:
			data = zoom(data, zoom_factor, order=1, mode='wrap')

		data -= np.median(data) # Zero offset
		data /= (1.4826*mad(data)) # Renormalize
	
		# Select region of pca0 associated with window position
		if self.wind_mode == 'WINDOW':
			x1 = self.x0; y1 = self.y0
		elif self.wind_mode == 'STRIPE':
			x1 = 0; y1 = self.y0
		else:
			x1 = 0; y1 = 0
	
		# print(y1, self.naxis2) This appears to be a stub
		x2 = x1 + self.naxis1
		y2 = y1 + self.naxis2
		# Make sure x2 and y2 are valid
		if (x2 > data.shape[0] or y2 > data.shape[1]):
			_log.warn('Specified window size does not fit within detector array!')
			_log.warn('X indices: [%s,%s]; Y indices: [%s,%s]; XY Size: [%s, %s]' % 
						(x1,x2,y1,y2,data.shape[0],data.shape[1]))
			os.sys.exit()
		self.pca0 = data[y1:y2,x1:x2]
		
		# How many reference pixels on each border?
		w = self.reference_pixel_border_width # Easier to work with
		lower = w-y1; upper = w-(det_size-y2)
		left = w-x1; right = w-(det_size-x2)
		ref_all = np.array([lower,upper,left,right])
		ref_all[ref_all<0] = 0
		self.ref_all = ref_all


	def message(self, message_text):
		"""
		Used for status reporting
		"""
		if self.verbose is True:
			print('NG: ' + message_text + ' at DATETIME = ', \
				  datetime.datetime.now().time())

        def interp(self, data):
                """
                Interpolation function to convert distribution into a callable
                function to extract read noise values from.

                Parameters:
                        data - Data to draw from
                """
                X = np.sort(data)
                x_range = np.linspace(0,1,len(X))
                func = UnivariateSpline(x_range, X, k=4, s=0)
                return func

        def make_det_vals(self, data, N):
                """
                Draw random samples from distribution to create read
                noise values.

                Parameters:
                        data - Data to draw from
                        N - number of times to sample
                """
                np.random.seed(1)
                rands = np.random.random(N)
                rands_sort = np.sort(rands)
                det_vals = self.interp(data)(rands_sort)
                np.random.seed(1)
                np.random.shuffle(det_vals)
                return det_vals      

	def white_noise(self, nstep=None):
		"""
		Generate white noise for an HxRG including all time steps
		(actual pixels and overheads).

		Parameters:
			nstep - Length of vector returned
		"""
		return(np.random.standard_normal(nstep))    

	def pink_noise(self, mode):
		"""
		Generate a vector of non-periodic pink noise.
  
		Parameters:
			mode - Selected from {'pink', 'acn'}
		"""

		# Configure depending on mode setting
		if mode is 'pink':
			nstep  = 2*self.nstep
			nstep2 = 2*self.nstep2 # JML
			f = self.f2
			p_filter = self.p_filter2
		else:
			nstep  = self.nstep
			nstep2 = self.nstep2 # JML
			f = self.f1
			p_filter = self.p_filter1

		# Generate seed noise
		mynoise = self.white_noise(nstep2)
  
		# Save the mean and standard deviation of the first
		# half. These are restored later. We do not subtract the mean
		# here. This happens when we multiply the FFT by the pinkening
		# filter which has no power at f=0.
		the_mean = np.mean(mynoise[:nstep2//2])
		the_std = np.std(mynoise[:nstep2//2])
  
		# Apply the pinkening filter.
		thefft = np.fft.rfft(mynoise)
		thefft = np.multiply(thefft, p_filter)
		result = np.fft.irfft(thefft)
		result = result[:nstep//2] # Keep 1st half of nstep

		# Restore the mean and standard deviation
		result *= the_std / np.std(result)
		result = result - np.mean(result) + the_mean
		  
		# Done
		return(result)



	def mknoise(self, o_file, dit=None,rd_noise=None, pedestal=None, c_pink=None,
				u_pink=None, acn=None, pca0_amp=None,
				reference_pixel_noise_ratio=None, ktc_noise=None,
				bias_offset=None, bias_amp=None):
		"""
		Generate a FITS cube containing only noise.

		Parameters:
			o_file   - Output filename
			dit      - tiem of single exposure
			pedestal - Magnitude of pedestal drift in electrons
			rd_noise - Standard deviation of read noise in electrons
			c_pink   - Standard deviation of correlated pink noise in electrons
			u_pink   - Standard deviation of uncorrelated pink noise in
					   electrons
			acn      - Standard deviation of alterating column noise in
					   electrons
			pca0     - Standard deviation of pca0 in electrons
			reference_pixel_noise_ratio - Ratio of the standard deviation of
										  the reference pixels to the regular
										  pixels. Reference pixels are usually
										  a little lower noise.                                          
			ktc_noise   - kTC noise in electrons. Set this equal to
						  sqrt(k*T*C_pixel)/q_e, where k is Boltzmann's
						  constant, T is detector temperature, and C_pixel is
						  pixel capacitance. For an H2RG, the pixel capacitance
						  is typically about 40 fF.
			bias_offset - On average, integrations start here in electrons. Set
						  this so that all pixels are in range.
			bias_amp    - A multiplicative factor that we multiply PCA-zero by
						  to simulate a bias pattern. This is completely
						  independent from adding in "picture frame" noise.

		Note1:
		Because of the noise correlations, there is no simple way to
		predict the noise of the simulated images. However, to a
		crude first approximation, these components add in
		quadrature.

		Note2:
		The units in the above are mostly "electrons". This follows convention
		in the astronomical community. From a physics perspective, holes are
		actually the physical entity that is collected in Teledyne's p-on-n
		(p-type implants in n-type bulk) HgCdTe architecture.
		"""

		self.message('Starting mknoise()')

		# ======================================================================
		#
		# NOISE PARAMETERS
		#
		# Read in noise from hsim config file.
		#
		# ======================================================================

                if dit > 120:
                        self.rd_noise  = config_data['systematics']['rd']       if \
                                rd_noise     is None else rd_noise
                else:
                        self.rd_noise = config_data['systematics']['rd_lowexp'] if \
                                rd_noise     is None else rd_noise
		self.pedestal  = config_data['systematics']['pedestal'] if \
                                pedestal     is None else pedestal
		self.c_pink    = config_data['systematics']['c_pink']   if \
                                c_pink       is None else c_pink
		self.u_pink    = config_data['systematics']['u_pink']   if \
                                u_pink       is None else u_pink
		self.acn       = config_data['systematics']['acn']      if \
                                acn          is None else acn
		self.pca0_amp  = config_data['systematics']['pca0_amp'] if \
                                pca0_amp     is None else pca0_amp

		# Change this only if you know that your detector is different from a
		# typical H2RG.
		self.reference_pixel_noise_ratio = config_data['systematics']['ref_ratio'] if \
			reference_pixel_noise_ratio is None else reference_pixel_noise_ratio

		# These are used only when generating cubes. They are
		# completely removed when the data are calibrated to
		# correlated double sampling or slope images. We include
		# them in here to make more realistic looking raw cubes.
		self.ktc_noise = config_data['systematics']['ktc_noise']        if \
                                ktc_noise   is None else ktc_noise 
		self.bias_offset = config_data['systematics']['bias_offset'].   if \
                                bias_offset is None else bias_offset
		self.bias_amp    = config_data['systematics']['bias_amp']       if \
                                bias_amp    is None else bias_amp

		# ======================================================================

		# Initialize the result cube. For up-the-ramp integrations,
		# we also add a bias pattern. Otherwise, we assume
		# that the aim was to simulate a two dimensional correlated
		# double sampling image or slope image.
		self.message('Initializing results cube')
		result = np.zeros((self.naxis3, self.naxis2, self.naxis1), \
						  dtype=np.float32)
		if self.naxis3 > 1:
			# Inject a bias pattern and kTC noise. If there are no reference pixels,
			# we know that we are dealing with a subarray. In this case, we do not
			# inject any bias pattern for now.
			#if self.reference_pixel_border_width > 0:
			#	bias_pattern = self.pca0*self.bias_amp + self.bias_offset
			#else:
			#	bias_pattern = self.bias_offset
			
			# Always inject bias pattern. Works for WINDOW and STRIPE (JML)
			bias_pattern = self.pca0*self.bias_amp + self.bias_offset
			
			# Add in some kTC noise. Since this should always come out
			# in calibration, we do not attempt to model it in detail.
			bias_pattern += self.ktc_noise * \
						 np.random.standard_normal((self.naxis2, self.naxis1))

			# Ensure that there are no negative pixel values. Data cubes
			# are converted to unsigned integer before writing.
			#bias_pattern = np.where(bias_pattern < 0, 0, bias_pattern)
			# Updated to conform to Python >=2.6. (JML)
			#bias_pattern[bias_pattern < 0] = 0
			# Actually, I think this makes the most sense to do at the very end (JML)
		
			# Add in the bias pattern
			for z in np.arange(self.naxis3):
				result[z,:,:] += bias_pattern


                # Make read noise from distribution. This is different for each pixel.
                if self.rd_noise > 0:
                        self.message('Generating rd_noise')
                        w = self.ref_all
                        r = self.reference_pixel_noise_ratio  # Easier to work with
                        
                        hdu = fits.open(self.rn_file)
                        rn_data = hdu[0].data * self.rd_noise

                        st = rn_data[-4:,:].flatten()
                        sl = rn_data[:,:4].flatten()
                        sr = rn_data[:,-4:].flatten()
                        sc = rn_data[4:-4,4:-4].flatten()
                        
                        for z in np.arange(self.naxis3):
                                here = np.zeros((self.naxis2, self.naxis1))

                                # Noisy reference pixels for each side of detector
                                if w[0] > 0: # lower
                                        here[:w[0],:] = r * self.make_det_vals(sc, w[0]*self.naxis1).reshape(\
                                                                (w[0],self.naxis1)) * np.random.standard_normal(\
                                                                (w[0],self.naxis1))
                                if w[1] > 0: # upper
                                        here[-w[1]:,:] = r * self.make_det_vals(sc, w[1]*self.naxis1).reshape(\
                                                                (w[1],self.naxis1)) * np.random.standard_normal(\
                                                                (w[1],self.naxis1))
                                if w[2] > 0: # left
                                        here[:,:w[2]] = r * self.make_det_vals(sc, w[2]*self.naxis2).reshape(\
                                                                (self.naxis2, w[2])) * np.random.standard_normal(\
                                                                (self.naxis2, w[2]))
                                if w[3] > 0: # right
                                        here[:,-w[3]:] = r * self.make_det_vals(sc, w[3]*self.naxis2).reshape(\
                                                                (self.naxis2, w[3])) * np.random.standard_normal(\
                                                                (self.naxis2, w[3]))

                                # Noisy regular pixels
                                if np.sum(w) > 0: # Ref. pixels exist in frame
                                        here[w[0]:self.naxis2-w[1],w[2]:self.naxis1-w[3]] = self.make_det_vals(sc, \
                                                                ((self.naxis2-w[0]-w[1])*(self.naxis1-w[2]-w[3]))).reshape(\
                                                                (self.naxis2-w[0]-w[1],self.naxis1-w[2]-w[3]))
                                else: # No Ref. pixels, so add only regular pixels
                                        here = self.make_det_vals(sc, self.naxis2*self.naxis1).reshape(\
                                                                (self.naxis2,self.naxis1))

                                # Add the noise in to the result
                                result[z,:,:] += here


		# Add correlated pink noise.
		if self.c_pink > 0:
			self.message('Adding c_pink noise')
			tt = self.c_pink * self.pink_noise('pink') # tt is a temp. variable
			tt = np.reshape(tt, (self.naxis3, self.naxis2+self.nfoh, \
								 self.xsize+self.nroh))[:,:self.naxis2,:self.xsize]
			for op in np.arange(self.n_out):
				x0 = op * self.xsize
				x1 = x0 + self.xsize
				# By default fast-scan readout direction is [-->,<--,-->,<--]
				# If reverse_scan_direction is True, then [<--,-->,<--,-->]
				# Would be nice to include option for all --> or all <--
				modnum = 1 if self.reverse_scan_direction else 0
				if np.mod(op,2) == modnum:
					result[:,:,x0:x1] += tt
				else:
					result[:,:,x0:x1] += tt[:,:,::-1]


			

		# Add uncorrelated pink noise. Because this pink noise is stationary and
		# different for each output, we don't need to flip it.
		if self.u_pink > 0:
			self.message('Adding u_pink noise')
			for op in np.arange(self.n_out):
				x0 = op * self.xsize
				x1 = x0 + self.xsize
				tt = self.u_pink * self.pink_noise('pink')
				tt = np.reshape(tt, (self.naxis3, self.naxis2+self.nfoh, \
								 self.xsize+self.nroh))[:,:self.naxis2,:self.xsize]
				result[:,:,x0:x1] += tt


		# Add ACN
		if self.acn > 0:
			self.message('Adding acn noise')
			for op in np.arange(self.n_out):

				# Generate new pink noise for each even and odd vector.
				# We give these the abstract names 'a' and 'b' so that we
				# can use a previously worked out formula to turn them
				# back into an image section.
				a = self.acn * self.pink_noise('acn')
				b = self.acn * self.pink_noise('acn')

				# Pick out just the real pixels (i.e. ignore the gaps)
				a = a[np.where(self.m_short == 1)]
				b = b[np.where(self.m_short == 1)]

				# Reformat into an image section. This uses the formula
				# mentioned above.
				acn_cube = np.reshape(np.transpose(np.vstack((a,b))),
									  (self.naxis3,self.naxis2,self.xsize))

				# Add in the ACN. Because pink noise is stationary, we can
				# ignore the readout directions. There is no need to flip
				# acn_cube before adding it in.
				x0 = op * self.xsize
				x1 = x0 + self.xsize
				result[:,:,x0:x1] += acn_cube


		# Add PCA-zero. The PCA-zero template is modulated by 1/f.
		if self.pca0_amp > 0:
			self.message('Adding PCA-zero "picture frame" noise')
			gamma = self.pink_noise(mode='pink')
			zoom_factor = self.naxis2 * self.naxis3 / np.size(gamma)
			gamma = zoom(gamma, zoom_factor, order=1, mode='mirror')
			gamma = np.reshape(gamma, (self.naxis3,self.naxis2))
			for z in np.arange(self.naxis3):
				for y in np.arange(self.naxis2):
					result[z,y,:] += self.pca0_amp*self.pca0[y,:]*gamma[z,y]
		


		# If the data cube has only 1 frame, reformat into a 2-dimensional
		# image.
		if self.naxis3 == 1:
			self.message('Reformatting cube into image')
			result = result[0,:,:]


		# If the data cube has more than one frame, convert to unsigned
		# integer
		if self.naxis3 > 1:
			# Data will be converted to 16-bit unsigned int
			# Ensure that there are no negative pixel values. 
			result[result < 0] = 0
			# And that anything higher than 65535 gets tacked to the top end
			result[result >= 2**16] = 2**16 - 1

			self.message('Converting to 16-bit unsigned integer')
			result = result.astype('uint16')
	   
		# Create HDU
		hdu = fits.PrimaryHDU(result)
		hdu.header.append()
		hdu.header.append(('RD_NOISE', self.rd_noise, 'Read noise'))
		hdu.header.append(('PEDESTAL', self.pedestal, 'Pedestal drifts'))
		hdu.header.append(('C_PINK', self.c_pink, 'Correlated pink'))
		hdu.header.append(('U_PINK', self.u_pink, 'Uncorrelated pink'))
		hdu.header.append(('ACN', self.acn, 'Alternating column noise'))
		hdu.header.append(('PCA0', self.pca0_amp, \
						   'PCA zero, AKA picture frame'))
		hdu.header['HISTORY'] = 'Created by NGHXRG version ' \
								+ str(self.nghxrg_version)
								
		# Write the result to a FITS file
		if o_file is not None:
			self.message('Writing FITS file')
			hdu.writeto(o_file, clobber='True')

		self.message('Exiting mknoise()')
		
		return hdu
