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

# Utils to format RST
def yield_sphinx_only_markup(lines):
    """
    Cleans-up Sphinx-only constructs (ie from README.rst),
    so that *PyPi* can format it properly.

    To check for remaining errors, install ``sphinx`` and run::

        python setup.py --long-description | sed -file 'this_file.sed' | rst2html.py  --halt=warning

    :param file_inp:     a `filename` or ``sys.stdin``?
    :param file_out:     a `filename` or ``sys.stdout`?`
    
    References
    ----------
    
    https://stackoverflow.com/questions/16367770/my-rst-readme-is-not-formatted-on-pypi-python-org

    Notes
    -----
    
    Check output with::
        
        python setup.py --long-description | rst2html.py > output.html

    """
    substs = [
        ## Selected Sphinx-only Roles.
        #
        (r':abbr:`([^`]+)`',        r'\1'),
        (r':ref:`([^`]+)`',         r'`\1`_'),
        (r':term:`([^`]+)`',        r'**\1**'),
        (r':dfn:`([^`]+)`',         r'**\1**'),
        (r':(samp|guilabel|menuselection):`([^`]+)`',        r'``\2``'),


        ## Sphinx-only roles:
        #        :foo:`bar`   --> foo(``bar``)
        #        :a:foo:`bar` XXX afoo(``bar``)
        #
        #(r'(:(\w+))?:(\w+):`([^`]*)`', r'\2\3(``\4``)'),
        (r':(\w+):`([^`]*)`', r'\1(``\2``)'),


        ## Sphinx-only Directives.
        #
        (r'\.\. doctest',           r'code-block'),
        (r'\.\. plot::',            r'.. '),
        (r'\.\. seealso',           r'info'),
        (r'\.\. glossary',          r'rubric'),
        (r'\.\. figure::',          r'.. '),


        ## Other
        #
        (r'\|version\|',              r'x.x.x'),
        
        ## added to make RADIS docs Pypi compatible
        (r'\.\. image::',          r'.. '),
#        (r'\.\. |CO2| replace:: CO\ :sub:`2`',          r'.. '),
    ]

    regex_subs = [ (re.compile(regex, re.IGNORECASE), sub) for (regex, sub) in substs ]

    def clean_line(line):
        try:
            for (regex, sub) in regex_subs:
                line = regex.sub(sub, line)
        except Exception as ex:
            print("ERROR: %s, (line(%s)"%(regex, sub))
            raise ex

        return line

    for line in lines:
        yield clean_line(line)

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
      packages=find_packages(),
      install_requires=[
				 'mpldatacursor',
#                        'numpy',          # let the user install it 
#                        'scipy',         # let the user install it  
#                        'matplotlib',    # let the user install it 
#                        'pandas',        # let the user install it 
#                        'sympy',         # let the user install it 
                        'pint>=0.7.2',  # Unit aware calculations
                        'publib>=0.2.1', # Plotting styles for Matplotlib
                        'plotly>=2.0.6',    # for line survey HTML output
                        'termcolor',     # terminal colors
                        'six',  # python 2-3 compatibility
                        'configparser', 
                        'astroquery',   # to fetch HITRAN databases
                        'json-tricks',   # to deal with non jsonable formats
                        'numpydoc',     # for Jedi (autocompletion) to recognize
                                        # numpy docstrings
                        'tables',       # for pandas to HDF5 export
                        'pytest',       # to run test suite
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