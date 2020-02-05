""" Install file for RADIS

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
 

"""
from __future__ import print_function
from __future__ import absolute_import
from setuptools import setup, find_packages
import codecs
import io
import re
import os
import io
from os.path import abspath, dirname, join, exists
import sys
from setuptips import yield_sphinx_only_markup

# Build description from README (PyPi compatible)
# (note: README.rst has been converted to README.md by register.py, and cleaned afterwards )
description = "A fast line-by-line code for high-resolution infrared molecular spectra"
readme_path = join(abspath(dirname(__file__)), "README.md")
if not exists(readme_path):
    long_description = description
else:
    with io.open(readme_path, encoding="utf-8") as f:
        long_description = f.read()

# Read version number from file
with open(join(dirname(__file__), "radis", "__version__.txt")) as version_file:
    __version__ = version_file.read().strip()

# Main install routine
setup(
    name="radis",
    version=__version__,
    description=description,
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/radis/radis",
    author="Erwan Pannier",
    author_email="erwan.pannier@gmail.com",
    license="GNU Lesser General Public License v3 (LGPLv3)",
    keywords=[
        "spectrum",
        "infrared",
        "spectra",
        "radiation",
        "nonequilibrium",
        "spectroscopy",
        "molecules",
        "HITRAN",
    ],
    packages=find_packages(),
    install_requires=[
        "numpy",
        "scipy",
        "matplotlib",
        "pandas",
        "plotly",
        "h5py",
        "numba",
        "mpldatacursor",
        "astropy",  # Unit aware calculations
        "pint>=0.7.2",  # Unit aware calculations
        "publib>=0.3.2",  # Plotting styles for Matplotlib
        "plotly>=2.5.1",  # for line survey HTML output
        "termcolor",  # terminal colors
        "six",  # python 2-3 compatibility
        "configparser",
        "astroquery>=0.3.9",  # to fetch HITRAN databases
        "json-tricks>=3.13.6",  # to deal with non jsonable formats
        "tables",  # for pandas to HDF5 export
        "pytest",  # to run test suite
        "h5py",  # to write HDF5 files
        "joblib",  # for parallel loading of SpecDatabase
        "numba",  # just-in-time compiler
    ],
    extras_require={
        "dev": [
            "numpydoc",  # for Jedi (autocompletion) to recognize
            "black",  # for code-linting in accordance to PEP8
            "isort",  # for sorting imports
        ]
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
        "Topic :: Scientific/Engineering",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Operating System :: OS Independent",
    ],
    include_package_data=True,  # add non .py data files in MANIFEST.in
    # package_data={'radis': ['radis/phys/units.txt']},
    zip_safe=False,  # impossible as long as we have external files read with __file__ syntax
    platforms="any",
)
