# -*- coding: utf-8 -*-
"""

@author: Erwan
"""

from __future__ import print_function, absolute_import, division, unicode_literals

import os
import sys
from os.path import dirname
import importlib

def getProjectRoot():
    ''' Return the full path of the project root'''

    return dirname(dirname(__file__))

def import_from_module(module, name):
    ''' Import object 'name' from module 'module' 
    raises AttributeError if name doesnt exist
    
    Parameters
    ----------
    
    module, name: str
        module path, object name
    
    '''
    impmodule = importlib.import_module(module)
    return getattr(impmodule , name)

class Chdir:
    ''' because we need to change directory to get into the RADIS folder to find
    the Git files and version, when imported from another program. This class
    then ensures we get back in the correct directory 
    
    Examples
    --------
    
    Do::
    
        cd = Chdir(os.path.dirname(__file__))
        try:
            (...)
        except:
            (...)
        finally:
            (...)
            cd.__del__()
    
    '''

    def __init__(self, newPath):
        self.savedPath = os.getcwd()
        os.chdir(newPath)

    def __del__(self):
        os.chdir(self.savedPath)

# %%
#==============================================================================
# Python 2/3 compatibility
#==============================================================================

if sys.version_info[0] == 2:
    import six.moves.configparser as configparser
    FileNotFoundError = OSError   # FileNotFoundError doesn't exist in Python 2....
    PermissionError = IOError   # PermissionError doesn't exist in Python 2....
else:
    import configparser
    FileNotFoundError = FileNotFoundError
    PermissionError = PermissionError
    
# %% Other stuff
    
class DatabankNotFound(FileNotFoundError):
    ''' Used when a line database is not found in radis.rc ''' 
    pass

# %%
#==============================================================================
# Optional packages
#==============================================================================

class NotInstalled(object):
    ''' A class to deal with optional packages 
    Will raise an error only if the package is used (but not if imported only)
    '''
    def __init__(self, name,info=''):
        self.__name = name
        self.__info = info

    def __getattr__(self,item):
        raise ImportError('The {0} package is required to use this '
                          'feature. {1}'.format(self.__name,self.__info))

    def __call__(self, *args, **kwargs):
        raise ImportError('The {0} package is required to use this '
                          'feature. {1}'.format(self.__name,self.__info))

# %% Test

def _test(*args, **kwargs):
    
    print('Project root:', getProjectRoot())
    
    return True

if __name__=='__main__':
    _test()
    
