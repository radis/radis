# -*- coding: utf-8 -*-
"""

Summary
-------

RADIS 

A code to simulate infrared spectra of molecules::

                *(((((((                                       
                 ((((((((((((              ,(((((              
                 ((((((((((((((((/   *((((((((((*              
                  ((((((((((((((((( ((((((((((((               
                      (((((((( (((((((((((((                   
                         *                                      
                       @@  *@@       ..  /@@                   
                  @@&  @@  *@@       @@  /@@  @@%              
              @@  @@&  @@  *@@  @@&  @@  /@@  @@%  @@          
              @@  @@&  @@  *@@  @@&  @@  /@@  @@%  @@          
              @@  @@&  @@  *@@  @@&  @@  /@@  @@%  @@  (@      
         ,@   @@  @@&  @@  *@@  @@&  @@  /@@  @@%  @@          
         @@   @@  @@&  @@  ,.                                  
                                    ,%&&&&&&&&&&&&&&&&&&&
          &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&      
           &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&       
             &&&&&&&&&&&&&&&&@@@@@@&@@@&&&@@@&&&&&&&&          
               &&&&&&&&&&&&&&&@@@@@@&&&&&&&&&&&&&&&            
                 &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&              
                   &&&&&&&&&&&&&&&&&&&&&&&&&&&.                
                       &&&&&&&&&&&&&&&&&&&                     
                               .**.                            
                                &&&,                         
                                 &&         

See Source code [1]_, Documentation [2]_, Package [3]_

Notes
-----

RADIS is nonequilibrium emission and absorption line-by-line code, for use by 
infrared spectroscopic that want to compare line databases, or experimentalist 
that want to fit their experimental line-of-sight spectra.

Written as a general purpose radiative solver, the code is built around the HITRAN, 
HITEMP and CDSD databases for molecules in their electronic ground state. Energy 
levels are read from tabulated databases or calculated from Dunham developments. 
Boltzmann, Treanor, and state specific vibrational distributions can be 
generated. A modular architecture makes it possible to add new species without 
modifications to the core code. Thus far, CO2, CO are featured for non-equilibrium 
calculations, and all species present in the HITRAN database are featured for 
equilibrium calculations. To fit experimental spectra, RADIS includes a line 
survey tool, an interface with a look-up database to improve fitting convergence 
times, and a multi-slab module with a radiative transfer equation solver to 
reproduce line-of-sight experiments. Validation cases against existing spectral 
codes and experimental results from various plasma sources are presented. 

The code will soon be available under under GNU General Public 
License v3.0 

References
----------

.. [1] Source code: `GitHub repository <https://github.com/radis/radis>`__

.. [2] Online Documentation: `Readthedocs.io <https://radis.readthedocs.io/en/latest/?badge=latest>`__

.. [3] Install as a package: `PyPi project <https://pypi.python.org/pypi/radis>`__

"""

from __future__ import absolute_import, division, print_function, unicode_literals

import os
from .misc.utils import getProjectRoot
from .misc.utils import Chdir as _chdir


# %% Global constants

DEBUG_MODE = False
'''bool: change this at runtime with::
    
    import radis
    radis.DEBUG_MODE = True

Use the printdbg() function in radis.misc, typically with::
    
    if __debug__: printdbg(...)

so that printdbg are removed by the Python preprocessor when running in
optimize mode::

    python -O *.py
'''

AUTO_UPDATE_SPEC = False
'''bool: experimental feature
used to autoupdate .spec files to the latest format, by simply saving
them again once they're loaded and fixed. 
Warning! Better have a copy of your files before that, or a way to regenerate 
them. 

Examples
--------

Add to the top of your script (once is enough!):: 
    
    import radis
    radis.AUTO_UPDATE_SPEC = True

See Also
--------
:func:`~radis.tools.database._update_to_latest_format`
'''

OLDEST_COMPATIBLE_VERSION = '1.0.0'
'''str: forces to regenerate cache files that were created in a previous version'''

# %% Version


def get_version(verbose=False, add_git_number=True):
    ''' Reads __version.txt__ and retrieve version number
    If ``add_git``, also appends Git commit number if we're on a gitted session '''

    # First get version
    with open(os.path.join(getProjectRoot(), '__version__.txt')) as version_file:
        version = version_file.read().strip()

    # Now get git info
    if add_git_number:
        import subprocess
        import sys
        cd = _chdir(os.path.dirname(__file__))
        try:
            label = subprocess.check_output('git describe')
            label = label.decode().strip()
            label = '-' + label
        except:
            if verbose:
                print("couldnt get git version: {0}".format(sys.exc_info()[1]))
            # probably not a git session. drop
            label = ''
        finally:
            version = version+label
            cd.__del__()

    return version


__version__ = get_version(add_git_number=False)


from .spectrum import *        # Spectrum object
from .io import *              # input / output
from .lbl import *             # line-by-line module
from .los import *             # line-of-sight module
from .phys import *            # conversion functions, blackbody objects
from .tools import *           # slit, database, line survey, etc.
from .test import *            # test
