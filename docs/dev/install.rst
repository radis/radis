.. _label_install:

Requirements
------------

RADIS has the following strict requirements:

- `python <https://www.python.org/>`_ 3
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
-------

You can either install RADIS from `pip`, the Python package manager. But if 
you want to modify the code and contribute, we suggest to clone the source 
from `GitHub <https://github.com/radis/radis>`_.  

**Use-only version** : cant modify the code

In a terminal, run::

    pip install --user radis

The 'pip' module has to be installed (by default if you've installed Python
with Anaconda). 

**Developer version**: to modify the code and contribute to the 
project. 

In a terminal, run::

    git clone https://github.com/radis/radis
    cd radis
    pip install -e .

The `-e` (editable) command creates a link from the local folder `./` folder into Python 
site-packages.

To make sure the install worked, run the :ref:`first example <label_first_example>`
from the Quick Start page. Then, you're all set. 


Test 
----

If you want to modify the source code, you need to ensure that you don't break
any of the existing tests. 
Refer to the :ref:`Test Section <label_dev_test>` to learn how to run the 
tests locally. 



Update 
------

With Pip you can keep the package up-to-date with::

    pip install radis --upgrade


In the developer version (installed with `pip -e`), use git to `pull` the latest changes from 
`GitHub <https://github.com/radis/radis>`_.


Help
----

If you encounter any problem, please `report an Issue <https://github.com/radis/radis/issues?utf8=%E2%9C%93&q=is%3Aissue>`_ on GitHub.  

You can also ask for advice on the `Q&A forum <https://groups.google.com/forum/#!forum/radis-radiation>`__ 
or the community chat:

.. image:: https://badges.gitter.im/Join%20Chat.svg
    :target: https://gitter.im/radis-radiation/community
    :alt: Gitter
