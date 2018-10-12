This module shows how to use the :class:`~radis.spectrum.spectrum.Spectrum` class, 
and the different methods that are associated: 
:py:meth:`~radis.spectrum.spectrum.Spectrum.rescale_path_length`,
:py:meth:`~radis.spectrum.spectrum.Spectrum.rescale_mole_fraction`, 
:py:meth:`~radis.spectrum.spectrum.Spectrum.apply_slit`, 
:py:meth:`~radis.spectrum.spectrum.Spectrum.store`, etc. 

.. toctree::
   :maxdepth: 2
   
   howto

---------------------------------------------------------------------



===========================
How to generate a Spectrum?
===========================

   
Calculate a Spectrum
--------------------

Usually a Spectrum object is the output from a line-by-line (LBL) radiation code. 
Refer to the LBL documentation for that. Example using the RADIS LBL module with
:py:func:`~radis.lbl.calc.calc_spectrum`::

    from radis import calc_spectrum 
    s = calc_spectrum(...)        # s is a Spectrum object 
    
Or with the :class:`~radis.lbl.factory.SpectrumFactory` class, which can 
be used to batch-generate multiple spectra using the same line database::

    from radis import SpectrumFactory
    sf = SpectrumFactory(...)
    s = sf.eq_spectrum(...)
    s2 = sf.non_eq_spectrum(...)

Note that the :class:`~radis.lbl.factory.SpectrumFactory` class has the 
:meth:`~radis.lbl.loader.DatabankLoader.init_database` method that 
automatically retrieves a spectrum from a 
:class:`~radis.tools.database.SpecDatabase` if you calculated it already, 
or calculates it and stores it if you didn't. Very useful for spectra that 
requires long computational times!

   
Initialize from text file or arrays 
-----------------------------------

Spectrum objects can also be generated from numpy arrays, or files. 

From numpy arrays, use :py:meth:`~radis.spectrum.spectrum.Spectrum.from_array` ::

    # w, T are two numpy arrays 
    from radis import Spectrum
    s = Spectrum.from_array(w, T, 'transmittance_noslit', 
                               waveunit='nm', unit='I/I0')
                               
              
From a file, use :py:meth:`~radis.spectrum.spectrum.Spectrum.from_txt` ::
                 
    # 'exp_spectrum.txt' contains a spectrum
    from radis import Spectrum
    s = Spectrum.from_txt('exp_spectrum.txt', 'radiance', 
                               waveunit='nm', unit='mW/cm2/sr/nm')

Convenience functions have been added to handle the usual cases: 
:func:`~radis.spectrum.spectrum.calculated_spectrum`, 
:func:`~radis.spectrum.spectrum.transmittance_spectrum` and
:func:`~radis.spectrum.spectrum.experimental_spectrum`::

    # w, T, I are numpy arrays for wavelength, transmittance and radiance
    from radis import calculated_spectrum, transmittance_spectrum, experimental_spectrum
    s1 = calculated_spectrum(w, I, wunit='nm', Iunit='W/cm2/sr/nm')     # creates 'radiance_noslit'  
    s2 = transmittance_spectrum(w, T, wunit='nm')                       # creates 'transmittance_noslit'
    s3 = experimental_spectrum(w, I, wunit='nm', Iunit='W/cm2/sr/nm')   # creates 'radiance'    
    
    
Generate a Blackbody (Planck) function object
---------------------------------------------

In RADIS you can either use the :func:`~radis.phys.blackbody.planck` and
:func:`~radis.phys.blackbody.planck_wn` functions that generate Planck
radiation arrays for wavelength and wavenumber, respectively. 

Or, you can use the :func:`~radis.phys.blackbody.sPlanck` function that
returns a :class:`~radis.spectrum.spectrum.Spectrum` object, with all 
the associated methods (add in a line-of-sight, compare, etc.)

Example::
    
    s = sPlanck(wavelength_min=3000, wavelength_max=50000,
                T=288, eps=1)
    s.plot()
   
   
Load from a .spec file 
----------------------

A ``.spec`` file contains all the Spectrum spectral quantities as well as the input 
conditions used to generate it. To retrieve it use the 
:func:`~radis.tools.database.load_spec` function::
    
    s = load_spec('my_spectrum.spec')

Sometimes the ``.spec`` file may have been generated under a compressed format 
where redundant spectral quantities have been removed (for instance, transmittance
if you already know absorbance). Use the :py:meth:`~radis.spectrum.spectrum.Spectrum.update` 
method to regenerate missing spectral quantities::

    s = load_spec('my_spectrum.spec', binary=True)
    s.update()    
    
If many spectra are stored in a folder, it may be time to set up a 
:class:`~radis.tools.database.SpecDatabase` structure to easily see all 
Spectrum conditions and get Spectrum that suits specific parameters 
 
    
    
====================================
How to access a Spectrum properties?
====================================

    
Get spectral quantities
-----------------------

Spectral quantities (see the list of quantities in :doc:`spectrum`) can be stored under 
different formats in a Spectrum object (with wavenumbers, or wavelengths
in air, or wavelengths in vacuum, for a given unit, etc.) 

It is recommended to use the :py:meth:`~radis.spectrum.spectrum.Spectrum.get` 
method to get exactly what you want::
    
    w, I = s.get('transmittance_noslit', wunit='cm-1')  
    _, R = s.get('radiance_noslit', wunit='nm', Iunit='W/cm2/sr/nm',
                 medium='air')  
        
The default quantities are::

    # Convoluted with slit function:
    'radiance', 'transmittance', 'emissivity'
    
    # Not convoluted: 
    'radiance_noslit', 'transmittance_noslit', 'emisscoeff', 'emisscoeff_continuum',
    'absorbance', 'abscoeff', 'abscoeff_continuum', 'emissivity_noslit'

See the latest list in the :data:`~radis.spectrum.utils.CONVOLUTED_QUANTITIES` 
and :data:`~radis.spectrum.utils.NON_CONVOLUTED_QUANTITIES`.
    
Get wavelength / wavenumber
---------------------------

Use the :py:meth:`~radis.spectrum.spectrum.Spectrum.get_wavelength` and
:py:meth:`~radis.spectrum.spectrum.Spectrum.get_wavenumber` methods::
    
    w_nm = s.get_wavelength()
    w_cm = s.get_wavenumber()

Print Spectrum conditions
-------------------------

Want to know under which calculation conditions was your Spectrum object 
generated, or under which experimental conditions it was measured? 
Just print it::

    print(s)
    
(that shows all spectral quantities stored in the object, all keys and 
values in the :attr:`~radis.spectrum.spectrum.Spectrum.conditions` dictionary, 
and all atoms/molecules stored in the :attr:`~radis.spectrum.spectrum.Spectrum.populations` 
dictionary)

You can also show the conditions only with 
:py:meth:`~radis.spectrum.spectrum.Spectrum.print_conditions`::

	s.print_conditions()

    

Plot spectral quantities
------------------------

Use :py:meth:`~radis.spectrum.spectrum.Spectrum.plot`::

    s.plot('radiance_noslit')
    
You can plot on the same figure as before using the convenient ``nfig`` parameter::
    
    s_exp.plot('radiance_noslit', nfig='same')

But for comparing different spectra you may want to use 
:py:func:`~radis.spectrum.compare.plot_diff` directly.
    
Plot the slit function that was applied
---------------------------------------

Use :py:meth:`~radis.spectrum.spectrum.Spectrum.plot_slit`. You can also
change the unit::
    
    s.apply_slit(0.5, 'cm-1')    # for instance
    s.plot_slit('nm')

    
Plot populations
----------------

Use :py:meth:`~radis.spectrum.spectrum.Spectrum.plot_populations`::

    s.plot_populations('vib', nunit='cm-3')
    

Save a Spectrum object
----------------------

To store use the :py:meth:`~radis.spectrum.spectrum.Spectrum.store` method::
    
    # s is a Spectrum object
    s.store('temp_file.spec')
    from radis import load_spec
    s2 = load_spec('temp_file.spec')
    assert s == s2  # compare both

The generated ``.spec`` file can be read (and edited) with any text editor. However, 
it may take a lot of space. If memory is important, you may use the ``compress=True`` 
argument which will remove redundant spectral quantities (for instance, transmittance
if you already know absorbance), and store the .spec file under binary format. Use
the :py:meth:`~radis.spectrum.spectrum.Spectrum.update` method to regenerate missing 
quantities::

    s.store('temp_file.spec', compress=True, if_exists_then='replace')
    s2 = load_spec('temp_file.spec')
    s2.update()    # regenerate missing quantities 

If calculating many spectra, you may want to manipulate all of them together 
in a :class:`~radis.tools.database.SpecDatabase`. You can add them to the existing 
database with the :py:meth:`~radis.tools.database.SpecDatabase.add` method::

    :class:`~radis.tools.database.SpecDatabase`
    
Note that if using the RADIS LBL code to generate your spectra, the 
:class:`~radis.lbl.factory.SpectrumFactory` class has the 
:meth:`~radis.lbl.loader.DatabankLoader.init_database` method that 
automatically retrieves a spectrum from a database if you calculated it already, 
or calculates it and stores it if you didn't. Very useful for spectra that 
requires long computational times!


Export to txt
-------------

Saving to .txt in general isn't recommended as you will loose some information (for instance,
the conditions). You better use :py:meth:`~radis.spectrum.spectrum.Spectrum.store` and export 
to .spec [hidden .json] format. 

If you really need to export a given spectral quantity to txt file (for use in another software, 
for instance), you can use the :py:meth:`~radis.spectrum.spectrum.Spectrum.savetxt` method that 
will export a given spectral quantity.

Example::

    s.savetxt('radiance_W_cm2_sr_um.csv', 'radiance_noslit', wunit='nm', Iunit='W/cm2/sr/µm')
    
    
=============================
How to manipulate a Spectrum?
=============================
   
    
Calculate missing quantities
----------------------------

Some spectral quantities can be infered from quantities stored in the Spectrum 
if enough conditions are given. For instance, transmittance can be recomputed
from the spectral absorption coefficient if the path length is stored in the 
conditions. 

The :py:meth:`~radis.spectrum.rescale.update` method can be used to do that. 
Example::

    # w, A are numpy arrays for wavenumber and absorption coefficient
    s = Spectrum.from_array(w, A, 'abscoeff', wunit='cm-1')
    s.update('transmittance_noslit')
    
Or, all derivable quantities can be computed using ``.update('all')`` or simply ``.update()``::

    s.update() 
    

Update Spectrum conditions
--------------------------

Spectrum conditions are stored in a :attr:`~radis.spectrum.spectrum.Spectrum.conditions` dictionary 

Conditions can be updated *a posteriori* by modifying the dictionary::

    s.conditions['path_length'] = 10    # cm 


Rescale Spectrum with new path length
-------------------------------------

Path length can be changed after the spectra was calculated with the 
:py:meth:`~radis.spectrum.spectrum.Spectrum.rescale_path_length` method. 
If the spectrum is not optically thin, this requires to solve the radiative 
transfer equation again, so the emisscoeff and abscoeff quantities 
will have to be stored in the Spectrum, or any equivalent combination 
(radiance_noslit and absorbance, for instance). 

Example:

    >>> from radis import load_spec
    >>> s = load_spec('co_calculation.spec')
    >>> s.rescale_path_length(0.5)      # calculate for new path_length
    
    
Rescale Spectrum with new mole fraction
---------------------------------------

.. warning::

    Rescaling mole fractions neglects the changes in collisional broadening

mole fraction can also be changed in post-processing, using the 
:py:meth:`~radis.spectrum.spectrum.Spectrum.rescale_mole_fraction` method
that works similarly to the :py:meth:`~radis.spectrum.spectrum.Spectrum.rescale_path_length` 
method. However, the broadening coefficients are left unchanged, which is 
valid for small mole fraction changes. However, for large mole fraction changes 
you will have to recalculate the spectrum from scratch.    
    
    >>> s.rescale_mole_fraction(0.02)   # calculate for new mole fraction


Apply instrumental slit function
--------------------------------

Use :py:meth:`~radis.spectrum.spectrum.Spectrum.apply_slit`::

    s.apply_slit(1.5)    # nm 
    
By default, convoluted spectra are thinner than non convoluted spectra, to remove 
side effects. Use the ``mode=`` argument to change this behaviour. 


    
Multiply, substract
-------------------

Sometimes you need to manipulate an experimental spectrum, to account
for calibration or remove a baseline. Spectrum operations are done 
just for that:

- :py:func:`~radis.spectrum.operations.add_constant`
- :py:func:`~radis.spectrum.operations.add_array`
- :py:func:`~radis.spectrum.operations.add_spectra`
- :py:func:`~radis.spectrum.operations.substract_spectra`

Most of these functions are implemented with the standard operators. Ex::
    
    ((s_exp - 0.1)*10).plot()   # works for a Spectrum s_exp

Note that these operators are purely algebraic and should not be used in place
of the line-of-sight functions, i.e, :py:func:`~radis.los.slabs.SerialSlabs` (``>``)
and :py:func:`~radis.los.slabs.MergeSlabs` (``//``)

Most of these functions will only work if there is only one 
`spectral quantity <https://radis.readthedocs.io/en/latest/spectrum/spectrum.html#spectral-quantities>`__
defined in the Spectrum. If there is any ambiguity, the following
functions can be used to discard all but one spectral quantity:
    
- :py:func:`~radis.spectrum.operations.Transmittance`
- :py:func:`~radis.spectrum.operations.Transmittance_noslit`
- :py:func:`~radis.spectrum.operations.Radiance`
- :py:func:`~radis.spectrum.operations.Radiance_noslit`

Offset, crop
------------

Use the associated functions: :py:func:`~radis.spectrum.operations.crop`,
:py:func:`~radis.spectrum.operations.offset`. 

They are also defined as methods of the Spectrum objects (see 
:py:meth:`~radis.spectrum.Spectrum.crop`, :py:meth:`~radis.spectrum.Spectrum.offset`),
so they can be used directly with::

    s.offset(3, 'nm') 
    s.crop(370, 380, 'nm') 

By default, using methods will modify the object inplace, using the functions will 
generate a new Spectrum. 
    

Remove a baseline
-----------------

Either use the :py:func:`~radis.spectrum.operations.add_constant` mentionned
above, which is implemented with the ``-`` operator::

    s2 = s - 0.1 
    
Or remove a linear baseline with:

- :py:func:`~radis.spectrum.operations.get_baseline`
- :py:func:`~radis.spectrum.operations.sub_baseline`


Calculate transmittance from radiance with Kirchoff's law
---------------------------------------------------------

RADIS can be used to infer spectral quantities from others if they can 
be derived. If on top that, equilibrium is assumed, then Kirchoff's law 
is used. See ``How to ... calculate missing quantities?`` and the 
:py:meth:`~radis.spectrum.rescale.update` method with argument ``assume_equilibrium=True``.
Example::

    s = calculated_spectrum(...)     # defines 'radiance_noslit')
    s.update('transmittance_noslit')
    s.plot('transmittance_noslit')
    
You can infer if a Spectrum is at (thermal) equilibrium with the 
:py:meth:`~radis.spectrum.spectrum.Spectrum.is_at_equilibrium` method, that 
looks up the declared spectrum conditions and ensures ``Tgas==Tvib==Trot``. 
It does not imply chemical equilibrium (mole fractions are still arbitrary)
    


===================================
How to manipulate multiple Spectra?
===================================

Build a line-of-sight profile
-----------------------------

RADIS allows the combination of Spectra such as::

    s_line_of_sight = (s_plasma_CO2 // s_plasma_CO) > (s_room_absorption) 

Refer to the `line-of-sight module <https://radis.readthedocs.io/en/latest/los/index.html>`__ 

.. _label_spectrum_howto_compare:
Compare two Spectra
-------------------

You can compare two Spectrum objects using the :py:meth:`~radis.spectrum.spectrum.Spectrum.compare_with` 
method, or simply the ``==`` statement (which is essentially the same thing)::

    s1 == s2 
    >>> True/False 
    s1.compare_with(s2)
    >>> True/False 
    
However, :py:meth:`~radis.spectrum.spectrum.Spectrum.compare_with` allows more freedom 
regarding what quantities to compare. ``==`` will compare everything of two spectra, 
including input conditions, units under which spectral quantities are stored, 
populations of species if they were saved, etc. In many situations we may want 
to simply compare the spectra themselves, or even a particular quantity like 
*transmittance_noslit*. Use::

    s1.compare_with(s2, spectra_only=True)                    # compares all spectral quantities 
    s1.compare_with(s2, spectra_only='transmittance_noslit')  # compares transmittance only 
    
The aforementionned methods will return a boolean array (True/False). If you 
need the difference, or ratio, or distance, between your two spectra, or simply 
want to plot the difference, you can use one of the predefined functions 
:func:`~radis.spectrum.compare.get_diff`, :func:`~radis.spectrum.compare.get_ratio`, 
:func:`~radis.spectrum.compare.get_distance`, :func:`~radis.spectrum.compare.get_residual` 
or the plot function :func:`~radis.spectrum.compare.plot_diff`::

    from radis.spectrum import plot_diff
    s1 = load_spec(temp_file_name)
    s2 = load_spec(temp_file_name2)
    plot_diff(s1, s2, 'radiance')

These functions usually require that the spectra are calculated on the same spectral 
range. When comparing, let's say, a calculated spectrum with experimental data, 
you may want to interpolate: you can have a look at the :py:meth:`~radis.spectrum.spectrum.Spectrum.resample` 
method. See :ref:`Interpolate a Spectrum on another <label_spectrum_howto_interpolate>` for details. 

In :func:`~radis.spectrum.compare.plot_diff`, you can choose to plot the absolute difference
(``method='diff'``), or the ratio (``method='ratio'``), or both::

    # Below we compare 2 CO2 spectra s_cdsd and s_hitemp previously calculated with two different line databases. 
    from radis import plot_diff
    plot_diff(s_cdsd, s_hitemp, method=['diff', 'ratio']) 

.. image:: ./cdsd4000_vs_hitemp_3409K.svg
    :alt: Figure CO2 CDSD-4000 vs HITEMP-2010



.. _label_spectrum_howto_interpolate:
Interpolate a Spectrum on another
---------------------------------

Let's interpolate a calculated spectrum on an experimental spectrum, 
using the :py:meth:`~radis.spectrum.spectrum.Spectrum.resample` and, for instance,
the :py:meth:`~radis.spectrum.spectrum.Spectrum.get_wavelength` method::

    # let's say we have two objects:
    s_exp = load_spec('...')
    s_calc = calc_spectrum(...)
    # resample:
    s_calc.resample(s_exp.get_wavelength(), 'nm')
    
Energy conservation is ensured and an error is raised if your interpolation 
is too bad. If you need to adjust the error threshold, see the parameters
in :py:meth:`~radis.spectrum.spectrum.Spectrum.resample`.



Create a database of Spectrum objects
-------------------------------------

Use the :class:`~radis.tools.database.SpecDatabase` class. It takes a 
folder as an argument, and converts it in a list of 
:class:`~radis.spectrum.spectrum.Spectrum` objects. Then, you can:

- see the properties of all spectra within this folder with
  :py:meth:`~radis.tools.database.SpecDatabase.see`
- select the Spectrum that match a given set of conditions with 
  :py:meth:`~radis.tools.database.SpecDatabase.get`, 
  :py:meth:`~radis.tools.database.SpecDatabase.get_unique` and 
  :py:meth:`~radis.tools.database.SpecDatabase.get_closest`
- fit an experimental spectrum against all precomputed spectra in 
  the folder with :py:meth:`~radis.tools.database.SpecDatabase.fit_spectrum`
    

    