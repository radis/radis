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
    pip install -e .[dev]

The `-e` (editable) command creates a link from the local folder `./` folder into Python 
site-packages.

To make sure the install worked, run the :ref:`first example <label_first_example>`
from the Quick Start page. Then, you're all set. 

Code linting
------------

Radis follows `Black <https://black.readthedocs.io/en/stable/>`__ style for code linting to
maintain consistent coding style across modules. Code style is checked using CI services
which run automatically on each pull request. **Black** is automatically installed when radis
is set-up in developer mode.

To format any file/files::

    black /path/to/file/or/directory/

You can include Black coding style `directly in some text editors <https://github.com/psf/black#editor-integration>`__

Alternatively, Black coding style can be checked automatically before each commit. For that all you need to do is to run the following command once::

    cd radis
    pre-commit install

On each commit, format will be fixed if it was incorrect. All you need to do is to commit a second time. Exemple::

    $ git commit -am "test"
    black....................................................................Failed
    - hook id: black
    - files were modified by this hook

    reformatted [ALL FAILING FILES]
    All done!
    1 file reformatted.
    
    $ git commit -am "test"
    black....................................................................Passed
    [develop XXX] test
     1 file changed, 1 insertion(+)

Note that pre-commit will always require you to commit again after a test was failed, because `it's safer <https://github.com/pre-commit/pre-commit/issues/532>`__. If for any reason you want to skip formatting you can commit with the ``--no-verify`` `argument <https://git-scm.com/docs/git-commit>`__.  



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
