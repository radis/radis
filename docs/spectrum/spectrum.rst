.. _label_spectrum:

===================
The Spectrum object
===================

RADIS has powerful tools to post-process spectra created by
:ref:`the line-by-line module <label_line_by_line>` or by other spectral codes.

At the core of the post-processing is the :class:`~radis.spectrum.spectrum.Spectrum` class,
which features methods to:

- :ref:`generate Spectrum objects<label_howto_generate_spectrum>` from text files or python arrays.
- :ref:`rescale a spectrum<label_spectrum_rescale>` without redoing the line-by-line calculation.
- apply :ref:`instrumental slit functions<label_spectrum_apply_slit>`.
- :ref:`plot<label_spectrum_plot>` with one line and in whatever unit.
- :ref:`crop, offset<label_spectrum_offset_crop>` or :ref:`interpolate <label_spectrum_howto_interpolate>`.
- :ref:`remove baselines<label_spectrum_remove_baseline>`.
- :ref:`multiply or add constants<label_spectrum_algebra>` as simply as with ``s=10*s`` or ``s=s+0.2`` in Python.
- :ref:`store<label_spectrum_store>` an experimental or a calculated spectrum while retaining the metadata.
- :ref:`compare different spectra<label_spectrum_howto_compare>`.
- combine multiple spectra along the :ref:`line-of-sight<label_spectrum_line_of_sight>`.
- manipulate a folder of spectra easily with :ref:`spectrum Databases<label_spectrum_database>`.
- compute transmittance from absorbance, or whatever :ref:`missing spectral arrays<label_spectrum_missing_quantities>`.
- use the :ref:`line survey<label_spectrum_linesurvey>` tool to identify each line.


.. minigallery:: radis.Spectrum


Refer to the guide below for an exhaustive list of all features:

---------------------------------------------------------------------

.. toctree::
   :maxdepth: 3

   spectrum

For any other question you can use the `Q&A forum <https://groups.google.com/forum/#!forum/radis-radiation>`__,
the `GitHub issues <https://github.com/radis/radis/issues>`__ or the
community chats on `Gitter <https://gitter.im/radis-radiation/community>`__ or
`Slack <https://radis.github.io/slack-invite/>`__ . |badge_gitter| |badge_slack|



---------------------------------------------------------------------


.. _label_howto_generate_spectrum:


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
    sf.fetch_databank("hitemp", load_columns='noneq')  # or 'hitran', 'exomol', etc.
    s = sf.eq_spectrum(...)
    s2 = sf.non_eq_spectrum(...)

Note that the :class:`~radis.lbl.factory.SpectrumFactory` class has the
:meth:`~radis.lbl.loader.DatabankLoader.init_database` method that
automatically retrieves a spectrum from a
:class:`~radis.tools.database.SpecDatabase` if you calculated it already,
or calculates it and stores it if you didn't. Very useful for spectra that
require long computational times!


Initialize from Python arrays
-----------------------------

The standard way to build a Radis Spectrum is from a dictionary of
:py:mod:`numpy` arrays::

    # w, k, I are numpy arrays for wavenumbers, absorption coefficient, and radiance.
    from radis import Spectrum
    s = Spectrum({"wavenumber":w, "abscoeff":k, "radiance_noslit":I},
                 units={"radiance_noslit":"mW/cm2/sr/nm", "abscoeff":"cm-1"})
Or::

    s = Spectrum({"abscoeff":(w,k), "radiance_noslit":(w,I)},
                 wunit="cm-1"
                 units={"radiance_noslit":"mW/cm2/sr/nm", "abscoeff":"cm-1"})

You can also use the :py:meth:`~radis.spectrum.spectrum.Spectrum.from_array`
convenience function::

    # w, T are two numpy arrays
    from radis import Spectrum
    s = Spectrum.from_array(w, T, 'transmittance_noslit',
                               wunit='nm', unit='') # adimensioned

Dimensionned arrays can also be used directly ::

    import astropy.units as u
    w = np.linspace(200, 300) * u.nm
    I = np.random.rand(len(w)) * u.mW/u.cm**2/u.sr/u.nm
    s = Spectrum.from_array(w, I, 'radiance_noslit')

Other convenience functions have been added to handle the usual cases:
:func:`~radis.spectrum.models.calculated_spectrum`,
:func:`~radis.spectrum.models.transmittance_spectrum` and
:func:`~radis.spectrum.models.experimental_spectrum`::

    # w, T, I are numpy arrays for wavelength, transmittance and radiance
    from radis import calculated_spectrum, transmittance_spectrum, experimental_spectrum
    s1 = calculated_spectrum(w, I, wunit='nm', Iunit='W/cm2/sr/nm')     # creates 'radiance_noslit'
    s2 = transmittance_spectrum(w, T, wunit='nm')                       # creates 'transmittance_noslit'
    s3 = experimental_spectrum(w, I, wunit='nm', Iunit='W/cm2/sr/nm')   # creates 'radiance'


Initialize from Specutils
-------------------------

Use :py:meth:`~radis.spectrum.spectrum.Spectrum.from_specutils` to convert
from a ``specutils`` :py:class:`specutils.spectra.spectrum1d.Spectrum1D` object ::

    from radis import Spectrum
    Spectrum.from_specutils(spectrum)


Initialize from a text file
---------------------------

Spectrum objects can also be generated directly from a text file.

From a file, use :py:meth:`~radis.spectrum.spectrum.Spectrum.from_txt` ::

    # 'exp_spectrum.txt' contains a spectrum
    from radis import Spectrum
    s = Spectrum.from_txt('exp_spectrum.txt', 'radiance',
                               wunit='nm', unit='mW/cm2/sr/nm')

It is, however, recommended to use the RADIS ``.spec`` json format to store
and load arrays :


Load from a .spec file
----------------------

A ``.spec`` file contains all the Spectrum spectral arrays as well as the input
conditions used to generate it. To retrieve it use the
:func:`~radis.tools.database.load_spec` function::

    s = load_spec('my_spectrum.spec')

Sometimes the ``.spec`` file may have been generated under a compressed format
where redundant spectral arrays have been removed (for instance, transmittance
if you already know absorbance). Use the :py:meth:`~radis.spectrum.spectrum.Spectrum.update`
method to regenerate missing spectral arrays::

    s = load_spec('my_spectrum.spec', binary=True)
    s.update()

If many spectra are stored in a folder, it may be time to set up a
:class:`~radis.tools.database.SpecDatabase` structure to easily see all
Spectrum conditions and get Spectrum that suits specific parameters

.. minigallery:: radis.load_spec


Load from a HDF5 file
---------------------

This is the fastest way to read a Spectrum object from disk. It keeps metadata
and units, and you can also load only a part of a very large spectrum.
Use :py:meth:`~radis.spectrum.spectrum.Spectrum.from_hdf5`  ::

    Spectrum.from_hdf5("rad_hdf.h5", wmin=2100, wmax=2200, columns=['abscoeff', 'emisscoeff'])




Calculate a test spectrum
-------------------------

You need a spectrum immediatly, to run some tests ? Use :py:func:`~radis.test.utils.test_spectrum` ::

    s = radis.test_spectrum()
    s.apply_slit(0.5, 'nm')
    s.plot('radiance')

This returns the CO spectrum from the :ref:`first documentation example <label_first_example>`

.. figure:: examples/co_spectrum_700K.png
:scale: 60 %


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

.. minigallery:: radis.sPlanck




.. _label_spectral_arrays:

Spectral Arrays
===============

A :class:`~radis.spectrum.spectrum.Spectrum` object can contain one spectral arrays, such as
``'radiance'`` for emission spectra, or ``'transmittance'`` for absorption spectra.
It can also contain both emission and absorption quantities to be later combined
with other spectra by solving the :ref:`radiative transfer equation <label_los_index>`.

Some variables represent quantities that have been convolved with an instrumental slit function,
as measured in experimental spectra:

- ``'radiance'``: the spectral radiance, convolved by the instrument function (typically in :math:``'mW/cm^2/sr/nm'``). This is sometimes confusingly called *spectral intensity*.
- ``'transmittance'``: the directional spectral transmittance (:math:``0`` to :math:``1``), convolved by the instrument function.
- ``'emissivity'``: the directional spectral emissivity (:math:``0`` to :math:``1``), convolved by the instrument function.
  The spectral emissivity is the radiance emitted by a surface
  divided by that emitted by a black body at the same temperature as that surface. This value is only
  defined under thermal equilibrium.

Other variables represent quantities that have not been convolved (theoretical spectra):

- ``'radiance_noslit'``: the spectral radiance (typically in :math:`mW/cm^2/sr/nm`). This
is sometimes confusingly called *spectral intensity*.
- ``'transmittance_noslit'``: the directional spectral transmittance (:math:``0`` to :math:``1``)
- ``'emissivity_noslit'``: spectral emissivity (``0`` to ``1``) *i.e.* the radiance emitted by a surface
  divided by that emitted by a black body at the same temperature as that surface. This value is only
  defined under thermal equilibrium.
- ``'emisscoeff'``: the directional spectral emission density (typically in :math:``'mW/cm^3/sr/nm'``).
- ``'absorbance'``: the directional spectral absorbance (no dimension), also called *optical depth*.
- ``'abscoeff'``: spectral absorption coefficient (typically in :math:``'cm^{-1}'``), also called *opacity*.
  This is sometimes found as the *extinction coefficient* in the literature (strictly speaking, *extinction*
  is *absorption* + *diffusion*, the latter being negligible in the infrared).
- ``'xsection'``: absorption cross-section, typically in cm2

Additionally, RADIS may calculate extra quantities such as:

- ``'emisscoeff_continuum'``: the pseudo-continuum part of the spectral emission density ``'emisscoeff'``, that can be
  generated by :class:`~radis.lbl.factory.SpectrumFactory`
- ``'abscoeff_continuum'`` the pseudo-continuum part of the spectral absorption coefficient ``'abscoeff'``, that can be
  generated by :class:`~radis.lbl.factory.SpectrumFactory`

See the latest list in the :data:`~radis.spectrum.utils.CONVOLUTED_QUANTITIES` and
:data:`~radis.spectrum.utils.NON_CONVOLUTED_QUANTITIES`.

Custom spectral arrays
----------------------

A :class:`~radis.spectrum.spectrum.Spectrum` object is built on top of a dictionary structure, and can handle
spectral arrays with any name.

Custom spectral arrays with arbitrary units can be defined when creating a Spectrum object, for instance::

    # w, I are two numpy arrays
    s = Spectrum.from_array(w, I, 'irradiance',
                               wunit='nm', unit='w/cm2/nm')

Although not recommended, it is also possible to directly edit the dictionary containing the objects.
For instance, this is done in
`CO2 radiative forcing example <https://github.com/radis/radis-examples/blob/master/ex_radiative_forcing_co2/radiative_forcing_co2.py>`__
to calculate irradiance from radiance (by multiplying by ``'pi'`` and changing the unit)::

    s._q['irradiance'] = s.get('radiance_noslit')[1]*pi
    s.units['irradiance'] = s.units['radiance_noslit'].replace('/sr', '')

The unit conversion methods will properly work with custom units.

.. warning::

    Rescaling or combining spectra with custom quantities may result in errors.

Relations between quantities
----------------------------

Most of the quantities above can be recomputed from one another. In a homogeneous
slab, one requires an emission spectral density, and an absorption spectral density,
to be able to recompute the other quantities (provided that conditions such as path length
are given). Under equilibrium, only one quantity is needed. Missing quantities
can be recomputed automatically with the :meth:`~radis.spectrum.spectrum.Spectrum.update`
method.

Units
-----

Default units are stored in the :attr:`~radis.spectrum.spectrum.Spectrum.units` dictionary

It is strongly advised not to modify the dictionary above. However, spectral arrays
can be retrieved in arbitrary units with the :meth:`~radis.spectrum.spectrum.Spectrum.get`
method.

When a spectral unit is convolved with :meth:`~radis.spectrum.spectrum.Spectrum.apply_slit`,
a new convolved spectral array is created. The unit of the convolved spectral array may be different,
depending on how the slit function was normalized. Several options are available in RADIS.
Please refer to the documentation of the :meth:`~radis.spectrum.spectrum.Spectrum.apply_slit` method.





How to access Spectrum properties?
==================================


Get spectral arrays
-------------------

Spectral Arrays of a Spectrum object can be stored in arbitrary
wavespace (wavenumbers, wavelengths in air, or wavelengths in vacuum) and arbitrary units.

Therefore, it is recommendeded to use the :py:meth:`~radis.spectrum.spectrum.Spectrum.get`
method to retrieve the quantity un the units you want::

    w, I = s.get('transmittance_noslit', wunit='cm-1')
    _, R = s.get('radiance_noslit', wunit='nm', Iunit='W/cm2/sr/nm',
                 medium='air')

Use with `return_units` to get dimensioned Astropy Quantities ::

    w, R  = s.get('radiance_noslit', return_units=True)
    # w, R are astropy quantities

See :ref:`spectral arrays <label_spectral _arrays>` for the list
of spectral arrays.

Get wavelength/wavenumber
-------------------------

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

(that shows all spectral arrays stored in the object, all keys and
values in the :attr:`~radis.spectrum.spectrum.Spectrum.conditions` dictionary,
and all atoms/molecules stored in the :attr:`~radis.spectrum.spectrum.Spectrum.populations`
dictionary)

You can also show the conditions only with
:py:meth:`~radis.spectrum.spectrum.Spectrum.print_conditions`::

	s.print_conditions()


.. _label_spectrum_plot:

Plotting
    -----

    The `plot` method can be used to visualize the resulting spectrum. Available plot types are:

    - 'absorption': plot absorption coefficient vs. wavenumber
    - 'transmittance': plot transmittance vs. wavenumber
    - 'radiance': plot radiance vs. wavenumber
    - 'intensity': plot spectral intensity vs. wavenumber
    - 'lines': plot individual spectral lines

    Each plot type also has additional optional parameters that can be passed to `plot_options`.

    Parameters
    ----------
    var : str, optional
        The type of plot to generate. Defaults to 'absorption'.
        Valid options are 'absorption', 'transmittance', 'radiance', 'intensity', and 'lines'.
    plot_options : dict, optional
        A dictionary of additional plot options. Valid keys and their descriptions are:

        - 'figsize' : (width, height) tuple in inches for the plot figure size. Defaults to (10, 6).
        - 'xlim' : (xmin, xmax) tuple in cm-1 for the x-axis limits. Defaults to (0, 5000).
        - 'ylim' : (ymin, ymax) tuple for the y-axis limits. Defaults to 'auto'.
        - 'linewidth' : float value for the line width. Defaults to 1.0.
        - 'color' : string or tuple value for the line color. Defaults to 'k'.
        - 'title' : string value for the plot title. Defaults to an automatically generated title.
        - 'xlabel' : string value for the x-axis label. Defaults to 'Wavenumber (cm$^{-1}$)'.
        - 'ylabel' : string value for the y-axis label. Defaults to 'Absorption Coefficient (cm$^{-1}$/(molecule cm$^{-2}$))' for 'absorption' plots,
                      'Transmittance' for 'transmittance' plots,
                      'Radiance (W/(cm$^{-1}$ sr))' for 'radiance' plots,
                      'Spectral Intensity (W/(cm$^{-1}$ sr))' for 'intensity' plots.

    def plot(self, var ='absorption', plot_options=None):

Plot spectral arrays
--------------------

Use :py:meth:`~radis.spectrum.spectrum.Spectrum.plot`::

    s.plot('radiance_noslit')

You can plot on the same figure as before using the convenient ``nfig`` parameter::

    s_exp.plot('radiance_noslit', nfig='same')

But for comparing different spectra you may want to use
:py:func:`~radis.spectrum.compare.plot_diff` directly.

Plot populations
----------------

Get or plot populations computed in calculations.
Use :py:meth:`~radis.spectrum.spectrum.Spectrum.get_populations`
or :py:meth:`~radis.spectrum.spectrum.Spectrum.plot_populations`::

    s.plot_populations('vib', nunit='cm-3')

.. minigallery:: radis.spectrum.spectrum.Spectrum.get_populations

.. _label_spectrum_linesurvey:

Plot line survey
----------------

Use the :py:meth:`~radis.spectrum.spectrum.Spectrum.line_survey` method.
Example of output::

    from radis import SpectrumFactory
    sf = SpectrumFactory(
                         wavenum_min=2380,
                         wavenum_max=2400,
                         mole_fraction=400e-6,
                         path_length=100,  # cm
                         isotope=[1],
                         )
    sf.load_databank('HITRAN-CO2-TEST')
    s = sf.eq_spectrum(Tgas=1500)
    s.apply_slit(0.5)
    s.line_survey(overlay='radiance_noslit', barwidth=0.01)


.. raw:: html

    <iframe id="igraph" src="//plotly.com/~erwanp/6.embed" width="650" height="420" seamless="seamless" scrolling="no"></iframe>


.. minigallery:: radis.spectrum.spectrum.Spectrum.line_survey



Know if a spectrum has nan
--------------------------

:py:meth:`~radis.spectrum.spectrum.Spectrum.has_nan` looks over all spectral
quantities. ``print(s)`` will also show the number of nans per quantity ::

    s = radis.test_spectrum()
    s.has_nan()


How to export ?
===============

.. _label_spectrum_store:

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
argument which will remove redundant spectral arrays (for instance, transmittance
if you already know absorbance), and store the .spec file under binary format. Use
the :py:meth:`~radis.spectrum.spectrum.Spectrum.update` method to regenerate missing
quantities::

    s.store('temp_file.spec', compress=True, if_exists_then='replace')
    s2 = load_spec('temp_file.spec')
    s2.update()    # regenerate missing quantities

If calculating many spectra, you may want to handle all of them together
in a :class:`~radis.tools.database.SpecDatabase`. You can add them to the existing
database with the :py:meth:`~radis.tools.database.SpecDatabase.add` method::

    db = SpecDatabase(r"path/to/database")     # create or loads database
    db.add(s)

Note that if using the RADIS LBL code to generate your spectra, the
:class:`~radis.lbl.factory.SpectrumFactory` class has the
:meth:`~radis.lbl.loader.DatabankLoader.init_database` method that
automatically retrieves a spectrum from a database if you calculated it already,
or calculates it and stores it if you didn't. Very useful for spectra that
requires long computational times!

Export to hdf5
--------------

This is the fastest way to dump a Spectrum object on disk (and also, it keeps metadata
and therefore units !). Use :py:meth:`~radis.spectrum.spectrum.Spectrum.to_hdf5`  ::

    s.to_hdf5('spec01.h5')


Export to txt
-------------

Saving to ``.txt`` in general isn't recommended as you will loose some information (for instance,
the conditions). You better use :py:meth:`~radis.spectrum.spectrum.Spectrum.store` and export
to ``.spec`` [a hidden ``.json``] format.

If you really need to export a given spectral arrays to txt file (for use in another software,
for instance), you can use the :py:meth:`~radis.spectrum.spectrum.Spectrum.savetxt` method that
will export a given spectral arrays::

    s.savetxt('radiance_W_cm2_sr_um.csv', 'radiance_noslit', wunit='nm', Iunit='W/cm2/sr/µm')

Export to Pandas
----------------

Use :py:meth:`~radis.spectrum.spectrum.Spectrum.to_pandas` ::

    s.to_pandas()

This will return a DataFrame with all spectral arrays as columns.


Export to Specutils
-------------------

Use :py:meth:`~radis.spectrum.spectrum.Spectrum.to_specutils` to convert
to a to ``specutils`` :py:class:`specutils.spectra.spectrum1d.Spectrum1D` object ::

    s.to_specutils()


.. minigallery:: radis.spectrum.spectrum.Spectrum.to_specutils



How to modify a Spectrum object?
================================

.. _label_spectrum_missing_quantities:

Calculate missing quantities
----------------------------

Some spectral arrays can be infered from quantities stored in the Spectrum
if enough conditions are given. For instance, transmittance can be recomputed
from the spectral absorption coefficient if the path length is stored in the
conditions.

The :py:meth:`~radis.spectrum.rescale.update` method can be used to do that.
In the example below, we recompute transmittance from the absorption coefficient
(opacity) ::

    # w, A are numpy arrays for wavenumber and absorption coefficient
    s = Spectrum.from_array(w, A, 'abscoeff', wunit='cm-1')
    s.update('transmittance_noslit')

All derivable quantities can be computed using ``.update('all')`` or simply ``.update()``::

    s.update()


Update Spectrum conditions
--------------------------

Spectrum conditions are stored in a :attr:`~radis.spectrum.spectrum.Spectrum.conditions` dictionary

Conditions can be updated *a posteriori* by modifying the dictionary::

    s.conditions['path_length'] = 10    # cm

.. _label_spectrum_rescale:

Rescale Spectrum with new path length
-------------------------------------

Path length can be changed after the spectra were calculated with the
:py:meth:`~radis.spectrum.spectrum.Spectrum.rescale_path_length` method.
If the spectrum is not optically thin, this requires solving the radiative
transfer equation again, so the ``emisscoeff`` and ``abscoeff`` (opacity) quantities
will have to be stored in the Spectrum, or any equivalent combination
(radiance_noslit and absorbance, for instance).

Example::

    from radis import load_spec
    s = load_spec('co_calculation.spec')
    s.rescale_path_length(0.5)      # calculate for new path_length


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


.. _label_spectrum_apply_slit:

Apply instrumental slit function
--------------------------------

Use :py:meth:`~radis.spectrum.spectrum.Spectrum.apply_slit`::

    s.apply_slit(1.5)    # nm

By default, convoluted spectra are thinner than non-convoluted spectra, to remove
side effects. Use the ``mode=`` argument to change this behavior.

.. minigallery:: radis.Spectrum.apply_slit


Plot the slit function that was applied
---------------------------------------

Use :py:meth:`~radis.spectrum.spectrum.Spectrum.plot_slit`. You can also
change the unit::

    s.apply_slit(0.5, 'cm-1')    # for instance
    s.plot_slit('nm')



.. _label_spectrum_algebra:

Multiply, subtract
------------------

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
`spectral arrays <https://radis.readthedocs.io/en/latest/spectrum/spectrum.html#spectral-quantities>`__
defined in the Spectrum. If there is any ambiguity, use the :py:meth:`~radis.spectrum.spectrum.Spectrum.take`
method. For instance, the following line is a valid RADIS command to plot
the spectral radiance of a spectrum with a low resolution::

    (10*(s.apply_slit(10, 'nm')).take('radiance')).plot()

Algebraic operations also work with dimensioned `~astropy.units.quantity.Quantity`.
For instance, remove a constant baseline in a given unit::

    s -= 0.1 * u.Unit('W/cm2/sr/nm')

The :py:meth:`~radis.spectrum.spectrum.Spectrum.max` function returns a dimensionned
value, therefore it can be used to normalize a spectrum directly :  ::

    s /= s.max()

Or below, we calibrate a Spectrum, assuming the spectrum units is in "count",
and that our calibration show we have 94 :math:`mW/cm2/sr/nm` per count.  ::

    s *= 94 * u.Unit("mW/cm2/sr/nm") / u.Unit("count")


.. _label_spectrum_offset_crop:

Offset, crop
------------

Use the associated functions: :py:func:`~radis.spectrum.operations.crop`,
:py:func:`~radis.spectrum.operations.offset`.

They are also defined as methods of the Spectrum objects (see
:py:meth:`~radis.spectrum.spectrum.Spectrum.crop`, :py:meth:`~radis.spectrum.spectrum.Spectrum.offset`),
so they can be used directly with::

    s.offset(3, 'nm')
    s.crop(370, 380, 'nm')

By default, using methods that will modify the object in place, using the functions will
generate a new Spectrum.


.. minigallery:: radis.Spectrum.crop


Normalize
---------

Use :py:meth:`~radis.spectrum.spectrum.Spectrum.normalize` directly, if your spectrum
only has one spectral arrays ::

    s.normalize()
    s.normalize(normalize_how="max")
    s.normalize(normalize_how="area")

You can also normalize only on limited range. Useful for noisy spectra ::

    s.normalize(wrange=(2250, 2500), wunit="cm-1", normalize_how="mean")

This returns a new spectrum and does not modify the Spectrum itself. To do so use::

    s.normalize(inplace=True)

.. _label_spectrum_chaining:

Chaining
--------

You can chain the various methods of Spectrum. For instance::

    s.normalize().plot()

Or::

    s.crop(4120, 4220, 'nm').apply_slit(3, 'nm').take('radiance')

If you want to create a new spectrum, don't forget to set ``inplace=False``
for the first command that allows it. i.e ::

    s2 = s.crop(4120, 4220, 'nm', inplace=False).apply_slit(3, 'nm').offset(1.5, 'nm')


.. _label_spectrum_remove_baseline:

Remove a baseline
-----------------

Either use the :py:func:`~radis.spectrum.operations.add_constant` mentionned
above, which is implemented with the ``-`` operator::

    s2 = s - 0.1

Or remove a linear baseline with:

- :py:func:`~radis.spectrum.operations.get_baseline`
- :py:func:`~radis.spectrum.operations.sub_baseline`

You could also use the functions available in :py:mod:`pyspecutils`,
see :py:meth:`~radis.spectrum.spectrum.Spectrum.to_specutils`.


Calculate transmittance from radiance with Kirchoff's law
---------------------------------------------------------

RADIS can be used to infer spectral arrays from others if they can
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


How to handle multiple Spectra?
===============================

.. _label_spectrum_line_of_sight:

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
populations of species if they were saved, etc. In many situations, we may want
to simply compare the spectra themselves, or even a particular quantity like
*transmittance_noslit*. Use::

    s1.compare_with(s2, spectra_only=True)                    # compares all spectral arrays
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

.. image:: cdsd4000_vs_hitemp_3409K.*
    :alt: https://radis.readthedocs.io/en/latest/_images/cdsd4000_vs_hitemp_3409K.svg


Plot in log scale
-----------------

If you wish to plot in a logscale, it can be done in the following way:'::

    fig, [ax0, ax1] = plot_diff(s_expe, s_test, normalize=False, verbose=False)
    ylim0 = ax0.get_ybound()
    ax0.set_yscale("log")
    ax0.set_ybound(ylim0)

Fit an experimental spectrum
----------------------------

RADIS does not include fitting algorithms. To fit an experimental spectrum, one
should use one of the widely available optimization algorithms from the Python
ecosystem, for instance :py:func:`scipy.optimize.minimize`.

The :py:func:`~radis.spectrum.compare.get_residual` and
:py:func:`~radis.spectrum.compare.get_residual_integral` functions can
be used to return a scalar to feed to the :py:func:`~scipy.optimize.minimize`
function.

A simple fitting procedure could be::

    from scipy.optimize import minimize
    from radis import calc_spectrum, experimental_spectrum

    s_exp = experimental_spectrum(...)

    def cost_function(T):
        calc_spectrum(Tgas=T,
                      ... # other parameters )
        return get_residual(s_exp, s)

    best = minimize(cost_function,
                      800, # initial value
                      bounds=[500, 2000],
                      )
    T_best = best.x

Note however that the performances of a fitting procedure can be
vastly improved by not reloading the line database every time. In that
case, it becomes interesting to use the :class:`~radis.lbl.factory.SpectrumFactory`
class.

An example of a script that uses the :class:`~radis.lbl.factory.SpectrumFactory`,
multiple fitting parameters, and plots the residual and the calculated spectrum
in real-time, can be found :ref:`in the Examples page <label_examples_multitemperature_fit>`


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



.. _label_spectrum_database:

Create a database of Spectrum objects
-------------------------------------

Use the :class:`~radis.tools.database.SpecDatabase` class. It takes a
folder as an argument, and converts it in a list of
:class:`~radis.spectrum.spectrum.Spectrum` objects. Then, you can:

- see the properties of all spectra within this folder with
  :py:meth:`~radis.tools.database.SpecList.see`
- select the Spectrum that match a given set of conditions with
  :py:meth:`~radis.tools.database.SpecList.get`,
  :py:meth:`~radis.tools.database.SpecList.get_unique` and
  :py:meth:`~radis.tools.database.SpecList.get_closest`
- fit an experimental spectrum against all precomputed spectra in
  the folder with :py:meth:`~radis.tools.database.SpecDatabase.fit_spectrum`

See more information about databases below.


.. include:: _database.rst




.. |badge_gitter| image:: https://badges.gitter.im/Join%20Chat.svg
                  :target: https://gitter.im/radis-radiation/community
                  :alt: Gitter

.. |badge_slack| image:: https://img.shields.io/badge/slack-join-green.svg?logo=slack
                  :target: https://radis.github.io/slack-invite/
                  :alt: Slack
