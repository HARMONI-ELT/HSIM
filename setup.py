#!/usr/bin/env python

'''setup.py code for HARMONI Simulator hsim

Written by Simon Zieleniewski

Last updated: 14-01-16

'''

from distutils.core import setup
import commands
import os
#Get SVN revision number using svnversion command
SVN = str(commands.getoutput("svnversion -c ./ | sed -e 's/[MS]//g' -e 's/^[[:digit:]]*://'"))


setup(name='HSIM',
      version=SVN,
      description='HARMONI Simulation Pipeline',
      author='Simon Zieleniewski & Sarah Kendrew',
      author_email='simon.zieleniewski@physics.ox.ac.uk',
      url='http://astroweb2.physics.ox.ac.uk/~zieleniewski/',
      py_modules = ['hsim'],
      packages=['src', 'src.modules'],
      data_files=[('Sim_data/Background/ASM_background/',
                   ['Sim_data/Background/ASM_background/radiance_resolution_0.15_angstroms_MICRONS.txt']),
                  ('Sim_data/Background/Detector/',
                   ['Sim_data/Background/Detector/KMOS_Detector_chip3_noise_prob_dist.txt']),
                  ('Sim_data/Throughput/ASM_throughput/',
                   ['Sim_data/Throughput/ASM_throughput/transmission_resolution_0.15_angstroms_MICRONS.txt']),
                  ('Sim_data/Throughput/',
                   ['Sim_data/Throughput/DetectorQE.txt',
                    'Sim_data/Throughput/EELT_mirror_reflectivity_mgf2agal.dat.txt',
                    'Sim_data/Throughput/HARMONIthruput.txt']),
                  ('Sim_data/R500/',
                   ['Sim_data/R500/HARMONI_R500_mode_data.txt']),
                  ('Sim_data/PSFs/LTAO/',
                   ['Sim_data/PSFs/LTAO/LTAOdata.txt']),
                  ('Sim_data/PSFs/SCAO/',
                   ['Sim_data/PSFs/SCAO/SCAOdata.txt']),
                  ('Sim_data/PSFs/GLAO/',
                   ['Sim_data/PSFs/GLAO/GLAOdata.txt']),
                  ('Output_cubes/', ['Output_cubes/info.txt'])]
      )
