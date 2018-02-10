from __future__ import print_function
from __future__ import absolute_import
from setuptools import setup, find_packages
import codecs
import io
import os
from os.path import dirname, join
import sys

# Check mandatory modules that we want users to install themselves
try:
    import numpy 
    import scipy
    import matplotlib
    import pandas
except ImportError:
    raise ImportError('Please install these librairies first (with Anaconda is '+\
                      'strongly recommended) \n >>> conda install numpy scipy '+\
                      'matplotlib pandas')

long_description = 'A non-equilibrium Radiative Solver for HITRAN-like database species '
if os.path.exists('README.rst'):
    long_description = codecs.open('README.rst', encoding="utf-8").read()

__version__ = '0.1.2'

setup(name='radis',
      version=__version__,
      description='A non-equilibrium Radiative Solver for HITRAN-like database species ',
    	long_description=long_description,
      url='https://github.com/radis/radis',
      author='Erwan Pannier',
      author_email='erwan.pannier@gmail.com',
      license='GNU Lesser General Public License v3 (LGPLv3)',
      packages=find_packages(),
      install_requires=[
				 'mpldatacursor',
#                        'numpy',          # let the user install it 
#                        'scipy',         # let the user install it  
#                        'matplotlib',    # let the user install it 
#                        'pandas',        # let the user install it 
#                        'sympy',         # let the user install it 
                        'pint>=0.7.2',  # Unit aware calculations
                        'publib>=0.1.11', # Plotting styles for Matplotlib
                        'plotly>=2.0.6',
                        'termcolor',     # terminal colors
                        'six',  # python 2-3 compatibility
                        'configparser', 
				],
      classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
        'Topic :: Scientific/Engineering',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.6',
        "Operating System :: OS Independent"],
	  include_package_data=True,
      #data_files = [('radis/phys', ['units.txt'])],
      zip_safe=True,
      platforms='any',
      # Plugins for special features 
      extras_require={
          'test': ['pytest']
          
      })