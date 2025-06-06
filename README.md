# HSIM3

HSIM3 is a dedicated pipeline for simulating observations with HARMONI on the Extremely Large Telescope. HSIM3 takes high spectral and spatial resolution input data cubes, encoding physical descriptions of astrophysical sources, and generates mock observed data cubes. The simulations incorporate detailed models of the sky, telescope, instrument, and detectors to produce realistic mock data ([Zieleniewski et al. 2015b](https://doi.org/10.1093/mnras/stv1860)).

HSIM3 is programmed in Python and the source code can be found at https://github.com/HARMONI-ELT/HSIM.

The purpose of HSIM is to perform in-depth simulations for several key science cases, as part of the design and development of the HARMONI integral field spectrograph. The goal is to be able to simulate the performance of the instrument folding in detailed parameters including the ELT AO performance, atmospheric effects and realistic detector statistics.

For updated information on how to run HSIM3, you can read the manual [hsim3.pdf](https://github.com/HARMONI-ELT/HSIM/blob/master/hsim/manual/hsim3.pdf)


## System Requirements
The pipeline is written in Python v3.10. We have tested it on Linux Ubuntu, although it should be possible to run on Windows and Mac OSX.

The required Python modules are:
- astropy 6.1.6
- numpy 2.1.3
- scipy 1.13.1
- matplotlib 3.9.2
- astro-tiptop 1.3.12

The code has been tested with the indicated package version, although more recent releases of these packages are likely to work as well.

## Tips & Tricks ##
We point out here several useful tips that have emerged from our extensive development process.

1. QFitsView is a great tool for visualising 3D datacubes: [Download Here](http://www.mpe.mpg.de/~ott/QFitsView/)

2. Datacubes can get very large very quickly and fill up the memory even on large computers. The input datacube can be tailored to suit the required needs. For example if one is only interested in the spatial information of an object (e.g. around an emission line) then the cube can be kept very small in the spectral dimension.

3. The flux density units for input datacubes require the usual flux (e.g., erg/s/cm^2/AA) to be divided by the spaxel scale in arcsec. This then gives e.g. erg/s/cm^2/AA/arcsec^2.


## Contact ##

For questions please contact miguel.pereira@iff.csic.es

Developers: Miguel Pereira Santaella, Laurence Routledge, Simon Zieleniewsk, Sarah Kendrew

https://elt.eso.org/instrument/HARMONI/

https://harmoni-elt.physics.ox.ac.uk/
