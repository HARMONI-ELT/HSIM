#!/usr/bin/env python

'''setup.py code for JWST/NIRSpec Simulator JSIM

Written by Simon Zieleniewski

Last updated: 01-12-15

'''

from distutils.core import setup
import commands
import os
#Get SVN revision number using svnversion command
SVN = str(commands.getoutput("svnversion -c ./ | sed -e 's/[MS]//g' -e 's/^[[:digit:]]*://'"))


setup(name='JSIM',
      version=SVN,
      description='JWST/NIRSpec IFU Simulation Pipeline',
      author='Simon Zieleniewski',
      author_email='simon.zieleniewski@physics.ox.ac.uk',
      url='http://astroweb2.physics.ox.ac.uk/~zieleniewski/',
      py_modules = ['jsim'],
      packages=['src', 'src.modules'],
      data_files=[('Sim_data/Background/',
                   ['Sim_data/Background/Zodiacal_background_cgs.csv']),
                  ('Sim_data/Throughput/',
                   ['Sim_data/Throughput/F070-G140M.csv',
                    'Sim_data/Throughput/F070-G140H.csv',
                    'Sim_data/Throughput/F100-G140M.csv',
                    'Sim_data/Throughput/F100-G140H.csv',
                    'Sim_data/Throughput/G235M.csv',
                    'Sim_data/Throughput/G235H.csv',
                    'Sim_data/Throughput/G395M.csv',
                    'Sim_data/Throughput/G395H.csv',
                    'Sim_data/Throughput/NIRSpecthruput.txt']),
                  ('Sim_data/Prism/',
                   ['Sim_data/Prism/Prism_res.csv']),
                  ('Output_cubes/', ['Output_cubes/info.txt'])]
      )
