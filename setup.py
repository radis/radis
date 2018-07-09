''' Install file for RADIS

Typical install procedure, plus:
    
- auto-convert README.rst to long_description, removing some sphinx-only syntax 
so it can be rendered by PyPi
- read version number from __version__.txt
- some dependencies should be installed manually (typically: all packages with 
compiled components such as numpy, pandas, etc.)


Examples
--------

Install (normal, use-only)::

    python setup.py install
    
Or (create an alias, so you can still edit)::
    
    python setup.py develop

Notes
-----

For developers:

when creating a new version, just update the __version__.txt file

to register it on Pypi see register.py::
    
    python register.py 
 

'''
from __future__ import print_function
from __future__ import absolute_import
from setuptools import setup, find_packages
import codecs
import io
import re
import os
from os.path import dirname, join
import sys
from setuptips import yield_sphinx_only_markup

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

# Build description from README (PyPi compatible)
description = 'A non-equilibrium Radiative Solver for HITRAN-like database species '
readme_lines = codecs.open('README.rst', encoding="utf-8").readlines()
long_description = ''.join(yield_sphinx_only_markup(readme_lines))

# Read version number from file
with open(join(dirname(__file__),'radis', '__version__.txt')) as version_file:
    __version__ = version_file.read().strip()
    
# Main install routine
setup(name='radis',
      version=__version__,
      description=description,
    	long_description=long_description,
      url='https://github.com/radis/radis',
      author='Erwan Pannier',
      author_email='erwan.pannier@gmail.com',
      license='GNU Lesser General Public License v3 (LGPLv3)',
      keywords=["spectrum", "infrared", "spectra", "radiation", "nonequilibrium"],
      packages=find_packages(),
      install_requires=[
				 'mpldatacursor',
#                        'numpy',          # let the user install it 
#                        'scipy',         # let the user install it  
#                        'matplotlib',    # let the user install it 
#                        'pandas',        # let the user install it 
#                        'sympy',         # let the user install it 
                        'pint>=0.7.2',  # Unit aware calculations
                        'publib>=0.3.1', # Plotting styles for Matplotlib
                        'plotly>=2.5.1',    # for line survey HTML output  
                        # TODO: update to >=2.5.1 as soon as https://github.com/plotly/plotly.py/issues/963 is fixed
                        'termcolor',     # terminal colors
                        'six',  # python 2-3 compatibility
                        'configparser', 
                        'astroquery',   # to fetch HITRAN databases
                        'json-tricks',   # to deal with non jsonable formats
                        'numpydoc',     # for Jedi (autocompletion) to recognize
                                        # numpy docstrings
                        'tables',       # for pandas to HDF5 export
                        'pytest',       # to run test suite
                        'h5py',         # to write HDF5 files
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
	  include_package_data=True,  # add non .py data files in MANIFEST.in
      package_data={'radis': ['radis/phys/units.txt']},
      zip_safe=False,  # impossible as long as we have external files read with __file__ syntax
      platforms='any')