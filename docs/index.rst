.. RADIS documentation master file, created by
   sphinx-quickstart on Tue Feb  6 03:32:15 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. |logo_png| image:: radis_ico.png

*****
RADIS
*****

A nonequilibrium infrared emission and absorption line-by-line code ... 
and a post-processing library to manipulate experimental and calculated spectra.


===============
Getting Started
===============

Install
-------

Assuming you have Python installed with the `Anaconda <https://www.anaconda.com/download/>`_ distribution just use::

    pip install radis 
    
**That's it!** You can now run your first example below.
If you encounter any issue, or to upgrade the package later, please refer to the 
:ref:`detailed installation procedure <label_install>` . 

.. _label_first_example:
Quick Start
-----------


Calculate a CO equilibrium spectrum from the HITRAN database, using the
:py:func:`~radis.lbl.calc.calc_spectrum` function. Output is a 
:ref:`Spectrum object <label_spectrum>`: ::

    from radis import calc_spectrum
    s = calc_spectrum(1900, 2300,         # cm-1
                      molecule='CO',
                      isotope='1,2,3',
                      pressure=1.01325,   # bar
                      Tgas=700,           # K
                      mole_fraction=0.1, 
                      path_length=1,      # cm
                      )
    s.apply_slit(0.5, 'nm')       # simulate an experimental slit
    s.plot('radiance')

.. figure:: examples/co_spectrum_700K.png
    :scale: 60 %

Calculate a CO *nonequilibrium* spectrum from the HITRAN database
(on your first call, this will calculate and cache the CO(X) rovibrational
energies): ::

    s2 = calc_spectrum(1900, 2300,         # cm-1
                      molecule='CO',
                      isotope='1,2,3',
                      pressure=1.01325,   # bar
                      Tvib=700,           # K
                      Trot=300,           # K
                      mole_fraction=0.1, 
                      path_length=1,      # cm
                      )
    s2.apply_slit(0.5, 'nm')
    s2.plot('radiance', nfig='same')    # compare with previous
    
The Quick Start examples automatically download the line databases from [HITRAN-2016]_, which is valid for temperatures below 700 K. 
For *high temperature* cases, you may need to use other line databases such as 
[HITEMP-2010]_ (typically T < 2000 K) or [CDSD-4000]_ (T < 5000 K). These databases must be described in a ``~/.radis`` 
:ref:`Configuration file <label_lbl_config_file>`. 

More complex :ref:`examples <label_examples>` will require to use the :py:class:`~radis.lbl.factory.SpectrumFactory`
class, which is the core of RADIS line-by-line calculations. 
:py:func:`~radis.lbl.calc.calc_spectrum` is a wrapper to :py:class:`~radis.lbl.factory.SpectrumFactory`
for the simple cases. 

Refer to the :ref:`Examples <label_examples>` section for more examples. 

In the browser
--------------

Alternatively, you can also run RADIS directly in the browser with the  
`RADIS Interactive Examples <https://github.com/radis/radis-examples#interactive-examples>`_ project. 
For instance, run the Quick Start example on the link below:

.. image:: https://mybinder.org/badge.svg 
    :target: https://mybinder.org/v2/gh/radis/radis-examples/master?filepath=first_example.ipynb
    :alt: https://mybinder.org/v2/gh/radis/radis-examples/master?filepath=first_example.ipynb

Or start a bare RADIS online session (no installation needed!):
    
.. image:: https://mybinder.org/badge.svg 
    :target: https://mybinder.org/v2/gh/radis/radis-examples/master?filepath=radis_online.ipynb
    :alt: https://mybinder.org/v2/gh/radis/radis-examples/master?filepath=radis_online.ipynb

   
---------------------------------------------------------------------

=======
Content
=======

.. toctree::
   :maxdepth: 2
   
   user
   developer
   references
   
   
Access Module Methods
---------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


