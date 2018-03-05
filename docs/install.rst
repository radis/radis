Requirements
============

RADIS has the following strict requirements:

- `python <https://www.python.org/>`_ 2.7 or 3.6
- `numpy <http://www.numpy.org/>`_
- `scipy <https://www.scipy.org/>`_ 
- `matplotlib <https://matplotlib.org/>`_
- `pandas <https://pandas.pydata.org/>`_ 

The above librairies require compiled components that usually fail to be 
installed with `pip`. It is enforced that RADIS users install them manually 
beforehand. We suggest using the `Anaconda <https://www.anaconda.com/download/>`_ 
distribution that contains all the required scientific packages above, plus 
many others. 

Other dependant packages will be installed by `pip` on the installation 
time:

- `mpldatacursor <https://github.com/joferkington/mpldatacursor>`_: used for plotting 
- `pint <https://pint.readthedocs.io>`_: for unit-aware calculations 
- `publib <https://github.com/erwanp/publib>`_: plotting styles for Matplotlib
- `plotly <https://plot.ly/>`_: used to produce html output graph
- `termcolor <https://pypi.python.org/pypi/termcolor>`_: because it's prettier
- `six <https://pypi.python.org/pypi/six>`_: for python 2-3 compatibility
- `configparser <https://pypi.python.org/pypi/configparser>`_: to read user config files
- `astroquery <https://astroquery.readthedocs.io/en/latest/>`_: to fetch HITRAN files


Install
=======

You can either install RADIS from `pip`, the Python package manager. But if 
you want to modify the code and contribute, we suggest to clone the source 
from Github.  

**Use-only version** (not recommended): cant modify the code

In a terminal, run::

    pip install --user radis

The 'pip' module has to be installed (by default if you've installed Python
with Anaconda). 

**Developer version** (recommended): to modify the code and contribute to the 
project. 

In a terminal, run::

    git clone https://github.com/radis/radis
    cd radis
    python setup.py develop --user

The `develop` command creates a link from your /radis folder into Python 
site-packages.


Test 
====

To make sure the install worked, run the following command from the console in
the ``radis\radis`` directory::

    pytest


Update 
======

With Pip you can keep the package up-to-date with::

    pip install radis --upgrade


In the developer version, use git to `pull` the latest changes from Github. 




