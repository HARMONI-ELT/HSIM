'''Misc utilities

Author: Simon Zieleniewski

Last updated: 23-06-15

'''

import time
import numpy as n
import astropy.io.fits as p
import sys



###Datafile path setup
##def path_setup(path):
##    '''function to setup paths to datafiles'''
##    import os
##    
##    script_path = os.path.abspath(__file__)
##    script_dir = script_path[:-len(__file__.split('/')[-1])]
##    out_path = os.path.join(script_dir, path)
##
##    return out_path

#Datafile path setup
def path_setup(path):
    '''function to setup paths to datafiles'''
    import os
    script_dir = os.path.split(os.path.abspath(__file__))[0]
    out_path = os.path.join(script_dir, os.path.relpath(path),'')

    return out_path



#Random number seed
def generate_seed(set_seed):
    '''function that generates random number seed used to
    generate noise estimates for datacubes.

    Inputs:

        set_seed: Accepted values:(0, 1-10). If 0, use time/date seed
                  If value between 1-10, use that value as seed.

    '''
    
    #Use current date (year, month, day) and time (hour, min, sec) for random number seed
    #If set_seed=True, use set seed value.
    set_seed = float(set_seed)
    if not set_seed:
        time_seed = time.gmtime()
        seednum = time_seed.tm_year+time_seed.tm_mon+time_seed.tm_mday+\
               time_seed.tm_hour+time_seed.tm_min+time_seed.tm_sec
        #n.random.seed(seednum)
        print 'SEEDNUM = ', seednum
        return seednum
    elif set_seed >= 1 and set_seed <= 10:
        #n.random.seed(int(set_seed))
        print 'SET SEED = ', int(set_seed)
        return int(set_seed)



#Update progress 
def update_progress(progress):
    barLength = 20 # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
##    if not isinstance(progress, float):
##        progress = 0
##        status = "error: progress var must be float\r\n"
##    if progress < 0:
##        progress = 0
##        status = "Halt...\r\n"
##    if progress >= 1:
##        progress = 1
##        status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = "\rPercent: [{0}] {1}% {2}".format( "="*block + " "*(barLength-block), progress*100, status)
    sys.stdout.write(text)
    sys.stdout.flush()



#Check if input cube spatial extent is smaller than
#one output spaxel
def point_source(cube, head, spax):
    '''Function that checks spatial extent
    of input cube compared to desired output.
    Pads out with zeros to equivalent size for
    a single spaxel (1D spectrum) if necessary.

    Inputs:

        - cube: input cube
        - head: header
        - spax: tuple containing spaxel scales [mas, mas]

    Outputs:

        - scube: scaled cube
        - head: updated header

    '''

    #Total extent in mas for x and y axes
    extent_x = head['NAXIS1']*head['CDELT1']
    extent_y = head['NAXIS2']*head['CDELT2']

    #check if either extent is less than 1 full output spaxel size in mas
    test_x = extent_x/float(spax[0])
    test_y = extent_y/float(spax[1])
#    if test_x < 1. or test_y < 1.:

    ratio_x = spax[0]/float(head['CDELT1'])
    ratio_y = spax[1]/float(head['CDELT2'])

    newsize_x = ratio_x
    newsize_y = ratio_y
    
    ncube = n.zeros((cube.shape[0],newsize_y,newsize_x),dtype=n.float32)
    ncube[:,ncube.shape[1]/2.-cube.shape[1]/2.+1:ncube.shape[1]/2.+cube.shape[1]/2.+1,
          ncube.shape[2]/2.-cube.shape[2]/2.+1:ncube.shape[2]/2.+cube.shape[2]/2.+1] = cube

    head['NAXIS1'] = int(newsize_x)
    head['NAXIS2'] = int(newsize_y)


    return ncube, head


#Check if input cube spatial extent is smaller than
#one output spaxel - pad to HARMONI FoV size
def fov_check(cube, head, spax, fov):
    '''Function that checks spatial extent
    of input cube compared to desired output.
    Pads out with zeros to equivalent size for
    HARMONI FoV if necessary.

    Inputs:

        - cube: input cube
        - head: header
        - spax: tuple containing spaxel scales [mas, mas]
        - fov: field of view (x, y) spaxels

    Outputs:

        - scube: scaled cube
        - head: updated header

    '''

    #Total extent in mas for x and y axes
    extent_x = head['NAXIS1']*head['CDELT1']
    extent_y = head['NAXIS2']*head['CDELT2']

    #check if either extent is less than 1 full output spaxel size in mas
    test_x = extent_x/float(spax[0])
    test_y = extent_y/float(spax[1])
##    if test_x < 1. or test_y < 1.:

    ratio_x = spax[0]/float(head['CDELT1'])
    ratio_y = spax[1]/float(head['CDELT2'])

    newsize_x = fov[0]*ratio_x
    newsize_y = fov[1]*ratio_y
    
    ncube = n.zeros((cube.shape[0],newsize_y,newsize_x),dtype=n.float32)
    ncube[:,ncube.shape[1]/2.-cube.shape[1]/2.+1:ncube.shape[1]/2.+cube.shape[1]/2.+1,
          ncube.shape[2]/2.-cube.shape[2]/2.+1:ncube.shape[2]/2.+cube.shape[2]/2.+1] = cube

    head['NAXIS1'] = int(newsize_x)
    head['NAXIS2'] = int(newsize_y)


    return ncube, head

