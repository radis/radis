# -*- coding: utf-8 -*-
"""
Summary
-------

:class:`~radis.spectrum.spectrum.Spectrum` class holder

Keeps all output data from a SpectrumFactory calculation. Allows to reapply a
different slit function in post-processing, and plot all spectral quantities
with any unit

Routine Listings
----------------

:func:`~radis.spectrum.spectrum.is_spectrum`


Examples
--------

Typical use::

    from radis calculated_spectrum
    s = calculated_spectrum(w, I, conditions={'case':'previously calculated by ##'})
    s.plot('radiance_noslit')
    s.apply_slit(0.5, shape='triangular')
    s.plot('radiance')

Spectrum objects can be modified, stored, resampled, rescaled, or retrieved after they have been created
:py:meth:`~radis.spectrum.spectrum.Spectrum.store`,
:py:meth:`~radis.spectrum.spectrum.Spectrum.rescale_path_length`,
:py:meth:`~radis.spectrum.spectrum.Spectrum.rescale_mole_fraction`,
:py:meth:`~radis.spectrum.spectrum.Spectrum.resample`,
:py:meth:`~radis.spectrum.spectrum.Spectrum.store`,
:py:func:`~radis.tools.database.load_spec` ::

    from radis import load_spec
    s = load_spec('co_calculation.spec')
    s.rescale_path_length(0.5)                  # calculate for new path_length
    s.rescale_mole_fraction(0.02)   # calculate for new mole fraction
    s.resample(w_new)               # resample on new wavespace
    s.store('co_calculation2.spec')

See the :ref:`Spectrum object <label_spectrum>` for more post-processing functions, or
how to generate a spectrum from text.

-------------------------------------------------------------------------------


"""

from copy import deepcopy
from os.path import basename
from warnings import warn

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from numpy import abs, diff
from publib import fix_style, set_style

# from radis.lbl.base import print_conditions
from radis.misc.arrays import count_nans, evenly_distributed, nantrapz
from radis.misc.debug import printdbg
from radis.misc.signal import resample
from radis.phys.air import air2vacuum, vacuum2air
from radis.phys.convert import cm2nm, conv2, nm2cm
from radis.phys.units import Unit, convert_universal
from radis.spectrum.rescale import rescale_mole_fraction, rescale_path_length, update
from radis.spectrum.utils import (
    CONVOLUTED_QUANTITIES,
    NON_CONVOLUTED_QUANTITIES,
    cast_waveunit,
    format_xlabel,
    make_up,
    make_up_unit,
    print_conditions,
    split_and_plot_by_parts,
)

# %% Spectrum class to hold results )


class Spectrum(object):
    """This class holds results calculated with the
    :py:class:`~radis.lbl.factory.SpectrumFactory` calculation, with other
    radiative codes, or experimental data. It can be used to plot different
    quantities a posteriori, or manipulate output units (for instance convert a
    spectral radiance per wavelength units to a spectral radiance per
    wavenumber).

    See more information on how to generate, edit or combine Spectrum objects
    on :ref:`the Spectrum object guide <label_spectrum>`.

    Parameters
    ----------
    quantities: dict of tuples   {'quantity':(wavenum, quantity)}
        where quantities are spectral quantities (absorbance, radiance, etc.)
        and wavenum is in cm-1
        example::

        >>> {'radiance_noslit:':(wavenum, radiance_noslit),
             'absorbance':(wavenum, absorbance)}
    units: dict
        units for quantities

    Other Parameters
    ----------------
    conditions: dict
        physical conditions and calculation parameters
    cond_units: dict
        units for conditions
    populations: dict
        a dictionary of all species, and levels. Should be compatible with other
        radiative codes such as Specair output. Suggested format:
        {molecules: {isotopes: {elec state: rovib levels}}}
        e.g::

            {'CO2':{1: 'X': df}}   # with df a Pandas Dataframe

    lines: pandas Dataframe
        all lines in databank (necessary for using
        :meth:`~radis.spectrum.spectrum.Spectrum.line_survey`). Warning if you want to
        play with the lines content: The signification of columns in `lines` may be
        specific to a database format. Plus, some additional columns may have been
        added by the calculation (e.g: `Ei` and `S` for emission integral and
        linestrength in SpectrumFactory). Refer to the code to know what they mean
        (and their units)
    wavespace: ``'nm'``, ``'cm-1'``, ``'nm_vac'`` or ``None``
        wavelength in air (``'nm'``), wavenumber (``'cm-1'``), or wavelength in vacuum (``'nm_vac'``).
        Quantities should be evenly distributed along this space for fast
        convolution with the slit function
        If ``None``, ``'wavespace'`` must be defined in ``conditions``.
        (non-uniform slit function is not implemented anyway... )
        Defaults None (but raises an error if wavespace is not defined in
        conditions neither)


    Other Parameters
    ----------------
    name: str, or None
        Give a name to this Spectrum object (helps debugging in multislab
        configurations). Default ``None``
    warnings: boolean
        if ``True``, test if inputs are valid, e.g, spectra are evenly distributed in
        wavelength, and raise a warning if not. Note that this take ~ 3.5 ms for
        a 20k points spectrum, when the rest of the creation process is only
        ~ 1.8ms (makes it 3 times longer, and can be a problem if hundreds of
        spectra are created in a row). Default ``True``


    Examples
    --------
    Manipulate a Spectrum calculated by RADIS::

        s = calc_spectrum(2125, 2300, Tgas=2000, databank='CDSD')
        s.print_conditions()
        s.plot('absorbance')
        s.line_survey(overlay='absorbance')
        s.plot('radiance_noslit', wunits='cm-1', Iunits='W/m2/sr/cm-1')
        s.apply_slit(5)
        s.plot('radiance')
        w, t = s.get('transmittance_noslit')  # for use in multi-slabs configs

    Any tuple of numpy arrays (w, I) can also be converted into a Spectrum object
    from the :class:`~radis.spectrum.spectrum.Spectrum` class directly, or using
    the :func:`~radis.spectrum.models.calculated_spectrum` function.
    All the following methods are equivalent::

        from radis import Spectrum, calculated_spectrum
        s1 = calculated_spectrum(w, I, wunit='nm', Iunit='mW/cm2/sr/nm')
        s2 = Spectrum.from_array(w, I, 'radiance_noslit',
                               waveunit='nm', unit='mW/cm2/sr/nm')
        s3 = Spectrum({'radiance_noslit': (w, I)},
                      units={'radiance_noslit':'mW/cm2/sr/nm'},
                      waveunit='nm')

    See more examples in the [Spectrum]_ page.

    Spectrum objects can be stored, retrieved, rescaled, resampled::

        from radis import load_spec
        s = load_spec('co_calculation.spec')
        s.rescale_path_length(0.5)                  # calculate for new path_length
        s.rescale_mole_fraction(0.02)   # calculate for new mole fraction
        s.resample(w_new)               # resample on new wavespace
        s.store('co_calculation2.spec')


    Notes
    -----
    Implementation:

        quantities are stored in ``self._q`` and ``self._q_conv`` dictionaries.
        They are better accessed with the :meth:`~radis.spectrum.spectrum.Spectrum.get`
        method that deals with units and wavespace

    Wavebase:

        quantites are stored either in wavenum or wavelength base, but this doesnt
        matter as they are retrieved / plotted with the
        :meth:`~radis.spectrum.spectrum.Spectrum.get` and :meth:`~radis.spectrum.spectrum.Spectrum.plot`
        methods which have units as input arguments


    Attributes
    ----------
    conditions : dict
        Stores computation / measurement conditions
    populations: dict
        Stores molecules, isotopes, electronic states and vibrational or
        rovibrational populations


    See Also
    --------

    :func:`~radis.spectrum.models.calculated_spectrum`,
    :func:`~radis.spectrum.models.transmittance_spectrum`,
    :func:`~radis.spectrum.models.experimental_spectrum`
    :meth:`~radis.spectrum.spectrum.Spectrum.from_array`
    :meth:`~radis.spectrum.spectrum.Spectrum.from_txt`
    :func:`~radis.tools.database.load_spec`


    References
    ----------

    .. [Spectrum] See the :ref:`Spectrum object page <label_spectrum>`
    """

    # hardcode attribute names, but can save a lot of memory if hundreds of spectra
    __slots__ = [
        "_q",
        "_q_conv",
        "units",
        "conditions",
        "cond_units",
        "populations",
        "lines",
        "name",
        "_slit",
        "file",
    ]

    def __init__(
        self,
        quantities,
        units=None,
        conditions=None,
        cond_units=None,
        populations=None,
        lines=None,
        waveunit=None,
        name=None,
        warnings=True,
    ):
        # TODO: make it possible to give quantities={'wavespace':w, 'radiance':I,
        # 'transmittance':T, etc.} directly too (instead of {'radiance':(w,I), etc.})

        # TODO: add help on creating a Spectrum from a dictionary

        # Check inputs
        # ... Replace None attributes with dictionaries
        if conditions is None:
            conditions = {}
        if units is None:
            units = {}
        if cond_units is None:
            cond_units = {}
        if populations is None:
            populations = {}
        self._init_annotations()  # typing hints for get() and plot()

        # Deal with deprecated inputs
        # ... wavespace renamed waveunit
        if "wavespace" in conditions:
            warn(
                DeprecationWarning("wavespace key in conditions is now named waveunit")
            )
            conditions["waveunit"] = conditions["wavespace"]
            del conditions["wavespace"]

        # Waveunit
        # ... Cast in standard waveunit format
        if waveunit is None and not "waveunit" in conditions:
            raise AssertionError(
                "waveunit ('nm', 'cm-1'?) has to be defined in `conditions`"
                + "or with waveunit="
            )
        if waveunit is not None:
            waveunit = cast_waveunit(waveunit)
        if "waveunit" in conditions:
            conditions["waveunit"] = cast_waveunit(conditions["waveunit"])
        # ... Make sure unit match
        if "waveunit" in conditions:
            if waveunit is not None and conditions["waveunit"] != waveunit:
                raise ValueError(
                    "waveunit defined in conditions ({0}) and directly ({1}) dont match".format(
                        conditions["waveunit"], waveunit
                    )
                )
        else:  # ... or define them in dictionary
            conditions["waveunit"] = waveunit

        # Check quantities format
        if len(quantities) == 0:
            raise AssertionError(
                "Spectrum is created with no quantities. Add "
                + "`radiance`, `transmittance` etc... with dict "
                + "format. e.g: {`radiance`: (w,I)}"
            )

        self._q = {}  #: dict: stores non convoluted quantities
        self._q_conv = {}  #: dict: stores convoluted (slit) quantities

        self._slit = {}  #: dict: hold slit function

        for k, v in quantities.items():
            try:
                assert len(v) == 2
            except AssertionError:
                raise AssertionError(
                    "Attributes should have format (wavenumber, "
                    + "quantity) . Error with `{0}`".format(k)
                )

            w, I = v

            # creates a copy of w,I
            self._add_quantity(k, w, I, warnings=warnings)

        # Finally, add our attributes
        self.conditions = conditions
        self.populations = populations
        self.lines = lines
        self.units = units
        """ dict: units for spectral quantities.
        """
        self.cond_units = cond_units
        self.name = name
        self.file = None  # used to store filename when loaded from a file

    # %% Constructors

    @classmethod
    def from_array(self, w, I, quantity, waveunit, unit, *args, **kwargs):
        """Construct Spectrum from 2 arrays.

        Parameters
        ----------

        w, I: array
            waverange and vector
        quantity: str
            spectral quantity name
        waveunit: ``'nm'``, ``'cm-1'``, ``'nm_vac'``
            unit of waverange:         wavelength in air (``'nm'``), wavenumber
            (``'cm-1'``), or wavelength in vacuum (``'nm_vac'``).
        unit: str
            spectral quantity unit (arbitrary). Ex: 'mW/cm2/sr/nm' for radiance_noslit
        *args, **kwargs
            see :class:`~radis.spectrum.spectrum.Spectrum` doc

        Other Parameters
        ----------------

        conditions: dict
            physical conditions and calculation parameters
        cond_units: dict
            units for conditions
        populations: dict
            a dictionary of all species, and levels. Should be compatible with other
            radiative codes such as Specair output. Suggested format:
            {molecules: {isotopes: {elec state: rovib levels}}}
            e.g::

                {'CO2':{1: 'X': df}}   # with df a Pandas Dataframe

        lines: pandas Dataframe
            all lines in databank (necessary for using
            :meth:`~radis.spectrum.spectrum.Spectrum.line_survey`). Warning if you want to
            play with the lines content: The signification of columns in `lines` may be
            specific to a database format. Plus, some additional columns may have been
            added by the calculation (e.g: `Ei` and `S` for emission integral and
            linestrength in SpectrumFactory). Refer to the code to know what they mean
            (and their units)

        Returns
        -------

        s: Spectrum
            creates a :class:`~radis.spectrum.spectrum.Spectrum` object

        Examples
        --------

        Create a spectrum::

            from radis import Spectrum
            s = Spectrum.from_array(w, I, 'radiance_noslit',
                                   waveunit='nm', unit='mW/cm2/sr/nm')

        To create a spectrum with absorption and emission components
        (e.g: ``radiance_noslit`` and ``transmittance_noslit``, or ``emisscoeff``
        and ``abscoeff``) call the :class:`~radis.spectrum.spectrum.Spectrum`
        class directly. Ex::

            from radis import Spectrum
            s = Spectrum({'abscoeff': (w, A), 'emisscoeff': (w, E)},
                         units={'abscoeff': 'cm-1', 'emisscoeff':'W/cm2/sr/nm'},
                         waveunit='nm')


        See Also
        --------

        :class:`~radis.spectrum.spectrum.Spectrum`,
        :func:`~radis.spectrum.models.calculated_spectrum`,
        :func:`~radis.spectrum.models.transmittance_spectrum`,
        :func:`~radis.spectrum.models.experimental_spectrum`
        :meth:`~radis.spectrum.spectrum.Spectrum.from_txt`
        :func:`~radis.tools.database.load_spec`,
        :ref:`the Spectrum page <label_spectrum>`
        """

        quantities = {quantity: (w, I)}
        units = {quantity: unit}

        # Update Spectrum conditions
        conditions = kwargs.pop("conditions", {})
        for k in ["path_length"]:  # usual conditions
            if k in kwargs:
                conditions[k] = kwargs.pop(k)

        return self(
            quantities, units, waveunit=waveunit, conditions=conditions, *args, **kwargs
        )

    @classmethod
    def from_txt(self, file, quantity, waveunit, unit, *args, **kwargs):
        """Construct Spectrum from txt file.

        Parameters
        ----------
        file: str
            file name
        quantity: str
            spectral quantity name
        waveunit: ``'nm'``, ``'cm-1'``, ``'nm_vac'``
            unit of waverange: wavelength in air (``'nm'``), wavenumber
            (``'cm-1'``), or wavelength in vacuum (``'nm_vac'``).
        unit: str
            spectral quantity unit
        *args, **kwargs
            the following inputs are forwarded to loadtxt: ``'delimiter'``, ``'skiprows'``
            The rest if forwarded to Spectrum. see :class:`~radis.spectrum.spectrum.Spectrum`
            doc


        Other Parameters
        ----------------
        delimiter: ``','``, etc.
            see :py:func:`numpy.loadtxt`
        skiprows: int
            see :py:func:`numpy.loadtxt`
        argsort: bool
            sorts the arrays in ``file`` by wavespace. Convenient way to load
            a file where points have been manually added at the end. Default ``False``.

        *Optional Spectrum parameters*

        conditions: dict
            physical conditions and calculation parameters
        cond_units: dict
            units for conditions
        populations: dict
            a dictionary of all species, and levels. Should be compatible with other
            radiative codes such as Specair output. Suggested format:
            {molecules: {isotopes: {elec state: rovib levels}}}
            e.g::

                {'CO2':{1: 'X': df}}   # with df a Pandas Dataframe

        lines: pandas Dataframe
            all lines in databank (necessary for using
            :meth:`~radis.spectrum.spectrum.Spectrum.line_survey`). Warning if you want to
            play with the lines content: The signification of columns in `lines` may be
            specific to a database format. Plus, some additional columns may have been
            added by the calculation (e.g: `Ei` and `S` for emission integral and
            linestrength in SpectrumFactory). Refer to the code to know what they mean
            (and their units)

        Returns
        -------
        s: Spectrum
            creates a :class:`~radis.spectrum.spectrum.Spectrum` object


        Examples
        --------
        Generate an experimental spectrum from txt. In that example the
        ``delimiter`` key is forwarded to :py:func:`~numpy.loadtxt`::

            from radis import Spectrum
            s = Spectrum.from_txt('spectrum.csv', 'radiance', waveunit='nm',
                                      unit='W/cm2/sr/nm', delimiter=',')


        To create a spectrum with absorption and emission components
        (e.g: ``radiance_noslit`` and ``transmittance_noslit``, or ``emisscoeff``
        and ``abscoeff``) call the :class:`~radis.spectrum.spectrum.Spectrum`
        class directly. Ex::

            from radis import Spectrum
            s = Spectrum({'abscoeff': (w, A), 'emisscoeff': (w, E)},
                         units={'abscoeff': 'cm-1', 'emisscoeff':'W/cm2/sr/nm'},
                         waveunit='nm')

        Notes
        -----
        Internally, the numpy :py:func:`~numpy.loadtxt` function is used and transposed::

            w, I = np.loadtxt(file).T

        You can use ``'delimiter'`` and '``skiprows'`` as arguments.

        See Also
        --------
        :func:`~radis.spectrum.models.calculated_spectrum`,
        :func:`~radis.spectrum.models.transmittance_spectrum`,
        :func:`~radis.spectrum.models.experimental_spectrum`,
        :meth:`~radis.spectrum.spectrum.Spectrum.from_array`,
        :func:`~radis.tools.database.load_spec`,
        :ref:`the Spectrum page <label_spectrum>`
        """

        # Get input for loadtxt
        kwloadtxt = {}
        for k in ["delimiter", "skiprows"]:
            if k in kwargs:
                kwloadtxt[k] = kwargs.pop(k)
        argsort = kwargs.pop("argsort", False)

        w, I = np.loadtxt(file, **kwloadtxt).T
        if argsort:
            b = np.argsort(w)
            w, I = w[b], I[b]
        quantities = {quantity: (w, I)}
        units = {quantity: unit}

        # Update Spectrum conditions
        conditions = kwargs.pop("conditions", {})
        for k in ["path_length"]:  # usual conditions
            if k in kwargs:
                conditions[k] = kwargs.pop(k)

        s = self(
            quantities, units, waveunit=waveunit, conditions=conditions, *args, **kwargs
        )

        # Store filename
        s.file = file
        return s

    # Public functions
    # %% ======================================================================
    # ----------------
    # XXX =====================================================================

    def get(self, var, wunit="default", Iunit="default", copy=True):
        """Retrieve a spectral quantity from a Spectrum object. You can select
        wavespace unit, intensity unit, or propagation medium.

        Parameters
        ----------
        var: variable ('absorbance', 'transmittance', etc.)
            Should be a defined quantity among :data:`~radis.spectrum.utils.CONVOLUTED_QUANTITIES`
            or :data:`~radis.spectrum.utils.NON_CONVOLUTED_QUANTITIES`.
            To get the full list of quantities defined in this Spectrum object use
            the :meth:`~radis.spectrum.spectrum.Spectrum.get_vars` method.
        wunit: ``'nm'``, ``'cm'``, ``'nm_vac'``.
            wavespace unit: wavelength in air (``'nm'``), wavenumber
            (``'cm-1'``), or wavelength in vacuum (``'nm_vac'``).
            if ``"default"``, default unit for waveunit is used. See
            :py:meth:`~radis.spectrum.spectrum.Spectrum.get_waveunit`.
        Iunit: unit for variable ``var``
            if ``"default"``, default unit for quantity `var` is used. See the
            :py:attr:`~radis.spectrum.spectrum.Spectrum.units` attribute.
            For ``var="radiance"``, one can use per wavelength (~ 'W/m2/sr/nm')
            or per wavenumber (~ 'W/m2/sr/cm-1') units

        Other Parameters
        ----------------
        copy: boolean
            if ``True``, returns a copy of the stored quantity (modifying it wont
            change the Spectrum object). Default ``True``.

        Returns
        -------
        w, I: array-like
            wavespace, quantity (ex: wavelength, radiance_noslit). For numpy
            users, note that these are copies (values) of the Spectrum quantity
            and not a view (reference): if you modify them the Spectrum is not
            changed

        Examples
        --------
        Get transmittance in cm-1::

            w, I = s.get('transmittance_noslit', wunit='cm-1')

        Get radiance (in wavelength in air)::

            _, R = s.get('radiance_noslit', wunit='nm', Iunit='W/cm2/sr/nm')

        See Also
        --------
        :meth:`~radis.spectrum.spectrum.Spectrum.get_radiance`,
        :meth:`~radis.spectrum.spectrum.Spectrum.get_radiance_noslit`,
        :ref:`the Spectrum page <label_spectrum>`
        """
        # TODO: allow to get radiance in (W/sr/cm2/nm) or (W/sr/cm2) multiplying
        # by FWHM.

        # check input
        if not var in self.get_vars():
            if var + "_noslit" in self.get_vars():
                raise KeyError(
                    ("`{0}` doesnt exist but `{1}` does".format(var, var + "_noslit"))
                    + ". Have you used .apply_slit()?"
                )
            else:
                raise KeyError(
                    "{0} not in quantity list: {1}".format(var, self.get_vars())
                )

        # Get quantity
        if var in CONVOLUTED_QUANTITIES:
            vartype = "convoluted"
            I = self._q_conv[var]
        elif var in NON_CONVOLUTED_QUANTITIES:  # non convoluted
            vartype = "non_convoluted"
            I = self._q[var]
        elif var in self._q:  # not declared, but exists. Assumes it's non_convoluted?
            vartype = "non_convoluted"
            I = self._q[var]
        else:
            raise ValueError("Unexpected quantity: {0}".format(var))

        # Copy
        if copy:
            I = I.copy()

        # Get wavespace (in correct unit, and correct medium)
        if wunit == "default":
            wunit = self.get_waveunit()
        else:
            wunit = cast_waveunit(wunit)
        if wunit == "cm-1":
            w = self.get_wavenumber(vartype, copy=copy)
        elif wunit == "nm":
            w = self.get_wavelength(medium="air", which=vartype, copy=copy)
        elif wunit == "nm_vac":
            w = self.get_wavelength(medium="vacuum", which=vartype, copy=copy)
        else:
            raise ValueError(wunit)

        # Convert y unit if necessary
        Iunit0 = self.units[var]
        if Iunit != "default" and Iunit != Iunit0:
            if var in ["radiance", "radiance_noslit"]:
                # deal with the case where we want to get a radiance in per
                # wavelength unit (~ W/sr/cm2/nm) in wavenumber units (~ W/sr/cm2/cm-1),
                # or the other way round
                I = convert_universal(
                    I,
                    Iunit0,
                    Iunit,
                    self,
                    per_nm_is_like="mW/sr/cm2/nm",
                    per_cm_is_like="mW/sr/cm2/cm-1",
                )
            elif var in ["emisscoeff"]:
                # idem for emisscoeff in (~ W/sr/cm3/nm) or (~ /sr/cm3/cm-1)
                I = convert_universal(
                    I,
                    Iunit0,
                    Iunit,
                    self,
                    per_nm_is_like="mW/sr/cm3/nm",
                    per_cm_is_like="mW/sr/cm3/cm-1",
                )
            elif var in ["absorbance"]:  # no unit
                assert Iunit in ["", "1"]
                # dont change the variable: I has no dimension
            elif var in ["transmittance"]:  # no unit
                assert Iunit in ["", "1"]
                # dont change the variable: I has no dimension
            else:
                I = conv2(I, Iunit0, Iunit)
        else:
            pass

        return w, I

    def _get_wavespace(self, which="any", copy=True):
        """Return wavespace (if the same for all quantities)

        Parameters
        ----------
        which: 'convoluted', 'non_convoluted', ``'any'``
            return wavelength for convoluted quantities, non convoluted quantities,
            or any. If ``any`` and both are defined, they have to be the same else
            an error is raised. Default ``any``.

        Other Parameters
        ----------------
        copy: boolean
            if ``True``, returns a copy of the stored waverange (modifying it wont
            change the Spectrum object). Default ``True``.


        Returns
        -------

        w: array_like
            (a copy of) spectrum wavespace for convoluted or non convoluted
            quantities
        """

        if which == "any":
            q_defined = "wavespace" in self._q
            q_conv_defined = "wavespace" in self._q_conv
            if q_defined and q_conv_defined:
                if not len(self._q["wavespace"]) == len(self._q_conv["wavespace"]):
                    raise ValueError(
                        "All wavespace not equal for calculated "
                        + "quantities. Can't use get_wavespace(). "
                        + "Specify which=`convoluted` or `non_convoluted`"
                    )
                if not np.allclose(self._q["wavespace"], self._q_conv["wavespace"]):
                    raise ValueError(
                        "All wavespace not equal for calculated "
                        + "quantities. Can't use get_wavespace(). "
                        + "Specify which=`convoluted` or `non_convoluted`"
                    )
                w = self._q["wavespace"]
            elif q_defined:
                w = self._q["wavespace"]
            else:
                w = self._q_conv["wavespace"]

        elif which == "convoluted":
            w = self._q_conv["wavespace"]

        elif which == "non_convoluted":
            w = self._q["wavespace"]

        else:
            raise ValueError(
                "which has to be one of 'convoluted', 'non_convoluted', "
                + "'any'. Got {0}".format(which)
            )

        if copy:
            w = w.copy()

        return w

    def get_wavelength(self, medium="air", which="any", copy=True):
        """Return wavelength in defined medium.

        Parameters
        ----------
        which: 'convoluted', 'non_convoluted', 'any'
            return wavelength for convoluted quantities, non convoluted quantities,
            or any. If any and both are defined, they have to be the same else
            an error is raised. Default any.
        medium: ``'air'``, ``'vacuum'``
            returns wavelength as seen in air, or vacuum. Default ``'air'``.
            See :func:`~radis.phys.air.vacuum2air`, :func:`~radis.phys.air.air2vacuum`

        Other Parameters
        ----------------
        copy: boolean
            if ``True``, returns a copy of the stored waverange (modifying it wont
            change the Spectrum object). Default ``True``.

        Returns
        -------
        w: array_like
            (a copy of) spectrum wavelength for convoluted or non convoluted
            quantities

        See Also
        --------
        :ref:`the Spectrum page <label_spectrum>`
        """

        # Check input
        if not medium in ["air", "vacuum"]:
            raise NotImplementedError("Unknown propagating medium: {0}".format(medium))

        # Now convert stored wavespace to the output unit
        w = self._get_wavespace(which=which, copy=copy)
        if self.get_waveunit() == "cm-1":
            w = cm2nm(w)  # vacuum wavelength

            # Correct for propagation medium (air, vacuum)
            if medium == "air":
                w = vacuum2air(w)
            else:  # medium == 'vacuum'
                pass  # no change needed

        elif self.get_waveunit() == "nm":  # nm air
            if medium == "air":
                pass  # no change needed
            else:
                w = air2vacuum(w)

        elif self.get_waveunit() == "nm_vac":  # nm vacuum
            if medium == "air":
                w = vacuum2air(w)
            else:
                pass  # no change needed

        else:
            raise ValueError(self.get_waveunit())

        return w

    def get_wavenumber(self, which="any", copy=True):
        """Return wavenumber (if the same for all quantities)

        Parameters
        ----------
        which: 'convoluted', 'non_convoluted', 'any'
            return wavenumber for convoluted quantities, non convoluted quantities,
            or any. If any and both are defined, they have to be the same else
            an error is raised. Default any.


        Other Parameters
        ----------------
        copy: boolean
            if ``True``, returns a copy of the stored waverange (modifying it wont
            change the Spectrum object). Default ``True``.

        Returns
        -------
        w: array_like
            (a copy of) spectrum wavenumber for convoluted or non convoluted
            quantities
        """
        w = self._get_wavespace(which=which, copy=copy)

        if self.get_waveunit() == "cm-1":  #
            pass

        elif self.get_waveunit() == "nm":  # wavelength air
            w = air2vacuum(w)
            w = nm2cm(w)

        elif self.get_waveunit() == "nm_vac":  # wavelength vacuum
            w = nm2cm(w)

        else:
            raise ValueError(self.get_waveunit())

        return w

    def get_radiance(self, Iunit="mW/cm2/sr/nm", copy=True):
        """Return radiance in whatever unit, and can even convert from ~1/nm to
        ~1/cm-1 (and the other way round)

        Other Parameters
        ----------------
        copy: boolean
            if ``True``, returns a copy of the stored waverange (modifying it wont
            change the Spectrum object). Default ``True``.


        See Also
        --------
        :meth:`~radis.spectrum.spectrum.Spectrum.get`,
        :meth:`~radis.spectrum.spectrum.Spectrum.get_radiance_noslit`,
        :ref:`the Spectrum page <label_spectrum>`
        """

        return self.get("radiance", Iunit=Iunit, copy=copy)[1]

    def get_radiance_noslit(self, Iunit="mW/cm2/sr/nm", copy=True):
        """Return radiance (non convoluted) in whatever unit, and can even
        convert from ~1/nm to ~1/cm-1 (and the other way round)

        Other Parameters
        ----------------
        copy: boolean
            if ``True``, returns a copy of the stored waverange (modifying it wont
            change the Spectrum object). Default ``True``.


        See Also
        --------
        :meth:`~radis.spectrum.spectrum.Spectrum.get`,
        :meth:`~radis.spectrum.spectrum.Spectrum.get_radiance`,
        :ref:`the Spectrum page <label_spectrum>`
        """

        return self.get("radiance_noslit", Iunit=Iunit, copy=copy)[1]

    def get_name(self):
        """Return Spectrum name.

        If not defined, returns either the
        :attr:`~radis.spectrum.spectrum.Spectrum.file` name if Spectrum was
        loaded from a file, or the ``'spectrum{id}'`` with
        the Python ``id`` object
        """
        try:
            self.name
        except AttributeError:
            warn(
                DeprecationWarning(
                    "Spectrum has no .name attribute and is probably outdated. Update!"
                )
            )
            self.name = None

        if self.name is not None:
            name = self.name
        elif self.file is not None:
            name = "{0}".format(basename(self.file))
        else:
            name = "spectrum{0}".format(id(self))

        return name

    def savetxt(self, filename, var, wunit="default", Iunit="default"):
        """Export spectral quantity var to filename.

        (note that by doing this you will loose additional information, such
         as the calculation conditions or the units. You better save a Spectrum
         object under a .spec file with :py:meth:`~radis.spectrum.spectrum.Spectrum.store`
         and load it afterwards with :py:func:`~radis.tools.database.load_spec`)

        Parameters
        ----------
        filename: str
            file name
        var: str
            which spectral variable ot export

        Other Parameters
        ----------------
        wunit, Iunit, medium: str
            see :meth:`~radis.spectrum.spectrum.Spectrum.get` for more information

        Notes
        -----
        Export variable as::

            np.savetxt(filename, np.vstack(self.get(var, wunit=wunit, Iunit=Iunit,
                                                medium=medium)).T, header=header)

        See Also
        --------
        :meth:`~radis.spectrum.spectrum.Spectrum.store`,
        :meth:`~radis.spectrum.spectrum.Spectrum.save`,
        :ref:`the Spectrum page <label_spectrum>`
        """

        # Get units to export
        if wunit == "default":
            wunit = self.get_waveunit()
        wunit = cast_waveunit(wunit)
        if wunit == "cm-1":
            xlabel = "Wavenumber (cm-1)"
        elif wunit == "nm":
            xlabel = "Wavelength [air] (nm)"
        elif wunit == "nm_vac":
            xlabel = "Wavelength [vacuum] (nm)"
        else:
            raise ValueError(wunit)

        if Iunit == "default":
            try:
                yunit = self.units[var]
            except KeyError:  # unit not defined in dictionary
                yunit = "a.u"
        else:
            yunit = Iunit

        header = "{0}\t{1} ({2})".format(xlabel, var, yunit)

        np.savetxt(
            filename,
            np.vstack(self.get(var, wunit=wunit, Iunit=Iunit)).T,
            header=header,
        )

    def update(self, quantity="all", optically_thin="default", verbose=True):
        """Calculate missing quantities: ex: if path_length and emisscoeff are
        given, recalculate radiance_noslit.

        Parameters
        ----------

        spec: Spectrum
        quantity: str
            name of the spectral quantity to recompute. If 'same', only the quantities
            in the Spectrum are recomputed. If 'all', then all quantities that can
            be derived are recomputed. Default 'all'.
        optically_thin: True, False, or 'default'
            determines whether to calculate radiance with or without self absorption.
            If 'default', the value is determined from the self_absorption key
            in Spectrum.conditions. If not given, False is taken. Default 'default'
            Also updates the self_absorption value in conditions (creates it if
            doesnt exist)

        See Also
        --------

        :ref:`the Spectrum page <label_spectrum>`
        """

        return update(
            self, quantity=quantity, optically_thin=optically_thin, verbose=verbose
        )

    # Rescale functions

    def rescale_path_length(
        self, new_path_length, old_path_length=None, inplace=True, force=False
    ):
        """Rescale spectrum to new path length. Starts from absorption
        coefficient and emission coefficient, and solves the RTE again for the
        new path length Convoluted values (with slit) are dropped in the
        process.

        Parameters
        ----------
        new_path_length: float
            new path length
        old_path_length: float, or None
            if None, current path length (conditions['path_length']) is used

        Other Parameters
        ----------------
        inplace: boolean
            if ``True``, modifies the Spectrum object directly. Else, returns
            a copy. Default ``True``.
        force: boolean
            if False, won't allow rescaling to 0 (not to loose information).
            Default ``False``

        Returns
        -------
        s: Spectrum
            Cropped Spectrum. If ``inplace=True``, Spectrum has been updated
            directly anyway.

        Notes
        -----
        Implementation:

            To deal with all the input cases, we first make a list of what has to
            be recomputed, and what has to be recalculated

        See Also
        --------
        :ref:`the Spectrum page <label_spectrum>`
        """

        return rescale_path_length(
            self,
            new_path_length=new_path_length,
            old_path_length=old_path_length,
            inplace=inplace,
            force=force,
        )

    def rescale_mole_fraction(
        self,
        new_mole_fraction,
        old_mole_fraction=None,
        inplace=True,
        ignore_warnings=False,
        force=False,
        verbose=True,
    ):
        """Update spectrum with new molar fraction Convoluted values (with
        slit) are dropped in the process.

        Parameters
        ----------
        new_mole_fraction: float
            new mole fraction
        old_mole_fraction: float, or None
            if None, current mole fraction (conditions['mole_fraction']) is used


        Other Parameters
        ----------------
        inplace: boolean
            if ``True``, modifies the Spectrum object directly. Else, returns
            a copy. Default ``True``.
        force: boolean
            if False, won't allow rescaling to 0 (not to loose information).
            Default ``False``

        Returns
        -------
        s: Spectrum
            Cropped Spectrum. If ``inplace=True``, Spectrum has been updated
            directly anyway.

        Notes
        -----
        Implementation:

            similar to rescale_path_length() but we have to scale abscoeff & emisscoeff
            Note that this is valid only for small changes in mole fractions. Then,
            the change in line broadening becomes significant

        See Also
        --------

        :ref:`the Spectrum page <label_spectrum>`
        """

        return rescale_mole_fraction(
            self,
            new_mole_fraction=new_mole_fraction,
            old_mole_fraction=old_mole_fraction,
            inplace=inplace,
            ignore_warnings=ignore_warnings,
            force=force,
            verbose=verbose,
        )

    def crop(self, wmin=None, wmax=None, wunit="default", inplace=True):
        """Crop spectrum to ``wmin-wmax`` range in ``wunit``   (inplace)

        Parameters
        ----------
        wmin, wmax: float, or None
            boundaries of spectral range (in ``wunit``)
        wunit: ``'nm'``, ``'cm-1'``, ``'nm_vac'``
            which waveunit to use for ``wmin, wmax``. If ``default``:
            use the default Spectrum wavespace defined with
            :meth:`~radis.spectrum.spectrum.Spectrum.get_waveunit`.


        Other Parameters
        ----------------
        inplace: bool
            if ``True``, modifies the Spectrum object directly. Else, returns
            a copy. Default ``True``.


        Returns
        -------
        s: Spectrum
            Cropped Spectrum. If ``inplace=True``, Spectrum has been updated
            directly anyway.


        Examples
        --------
        Crop to experimental Spectrum, and compare::

            from radis import calc_spectrum, load_spec, plot_diff
            s = calc_spectrum(...)
            s_exp = load_spec('typical_result.spec')
            s.crop(s_exp.get_wavelength.min(), s_exp.get_wavelength.max(), 'nm')
            plot_diff(s_exp, s)


        See Also
        --------
        :func:`radis.spectrum.operations.crop`,
        :func:`~radis.los.slabs.MergeSlabs`: if used with ``resample='full',
        out='transparent'``, this becomes the opposite of cropping: can be used
        to combine 2 adjacent spectra in one.
        """

        from radis.spectrum.operations import crop

        if wunit == "default":
            wunit = self.get_waveunit()

        return crop(self, wmin=wmin, wmax=wmax, wunit=wunit, inplace=inplace)

    def offset(self, offset, unit, inplace=True):
        # type: (Spectrum, float, str) -> Spectrum
        """Offset the spectrum by a wavelength or wavenumber  (inplace)

        Parameters
        ----------
        offset: float
            Constant to add to all quantities in the Spectrum.
        unit: 'nm' or 'cm-1'
            unit for ``offset``

        Other Parameters
        ----------------
        inplace: bool
            if ``True``, modifies the Spectrum object directly. Else, returns
            a copy. Default ``True``.

        Returns
        -------
        s: Spectrum
            Offset Spectrum. If ``inplace=True``, Spectrum has been updated
            directly anyway.

        See Also
        --------
        :func:`radis.spectrum.operations.offset`,
        :ref:`the Spectrum page <label_spectrum>`
        """

        from radis.spectrum.operations import offset as offset_func

        return offset_func(self, offset, unit, inplace=inplace)

    def get_integral(self, var, wunit="default", Iunit="default", **kwargs):
        """Returns integral of variable 'var' over waverange.

        Parameters
        ----------
        var: str
            spectral quantity to integate
        wunit: str
            over which waverange to integrated. If ``default``,
            use the default Spectrum wavespace defined with
            :meth:`~radis.spectrum.spectrum.Spectrum.get_waveunit`.
        Iunit: str
            default ``'default'``

            .. warning::
                this is the unit of the quantity, not the unit of the integral.
                Don't forget to multiply by ``wunit``

        Other Parameters
        ----------------
        kwargs: **dict
            forwarded to :meth:`~radis.spectrum.spectrum.Spectrum.get`

        Returns
        -------
        integral: float
            integral in [Iunit]*[wunit]

        See Also
        --------
        :meth:`~radis.spectrum.spectrum.Spectrum.get_power`,
        :ref:`the Spectrum page <label_spectrum>`
        """

        w, I = self.get(var, wunit=wunit, Iunit=Iunit, **kwargs)
        return abs(np.trapz(I, x=w))

    def get_power(self, unit="mW/cm2/sr"):
        """Returns integrated radiance (no slit) power density.

        Parameters
        ----------
        Iunit: str
            power unit.

        Returns
        -------
        P: float
            radiated power in ``unit``

        See Also
        --------
        :meth:`~radis.spectrum.spectrum.Spectrum.get_integral`,
        :ref:`the Spectrum page <label_spectrum>`
        """

        P = self.get_integral("radiance_noslit", wunit="nm", Iunit="mW/cm2/sr/nm")
        # P is in mW/cm2/sr/nm * nm
        return conv2(P, "mW/cm2/sr", unit)

    # %% Plotting routines

    def get_vars(self, which="any"):
        """Returns all spectral quantities stored in this object (convoluted or
        non convoluted)

        Parameters
        ----------
        which: 'any', 'convoluted', 'non convoluted'
        """
        if which == "any":
            varlist = list(self._q.keys()) + list(self._q_conv.keys())
        elif which == "convoluted":
            varlist = list(self._q_conv.keys())
        elif which == "non convoluted":
            varlist = list(self._q.keys())
        else:
            raise ValueError("Unexpected value: {0}".format(which))

        # remove wavespace
        varlist = [k for k in varlist if k != "wavespace"]
        return varlist

    def get_quantities(self, which="any"):
        """Returns all spectral quantities stored in this object (convoluted or
        non convoluted). Wrapper to
        :py:meth:`~radis.spectrum.spectrum.get_vars`

        Parameters
        ----------
        which: 'any', 'convoluted', 'non convoluted'
        """

        return self.get_vars(which=which)

    def _get_items(self):
        """Return a dictionary of tuples, e.g::

            {'radiance':(w,I), 'transmittance_noslit':(w_ns,T)}

        In the general case, users should use the
        :meth:`~radis.spectrum.spectrum.Spectrum.get` method

        .. warning::
            real quantities are returned, not copies.
        """
        items = {
            k: (self._q["wavespace"], v) for k, v in self._q.items() if k != "wavespace"
        }
        items.update(
            {
                k: (self._q_conv["wavespace"], v)
                for k, v in self._q_conv.items()
                if k != "wavespace"
            }
        )
        return items

    def plot(
        self,
        var=None,
        wunit="default",
        Iunit="default",
        show_points=False,
        nfig=None,
        yscale="linear",
        show_medium="vacuum_only",
        normalize=False,
        force=False,
        plot_by_parts=False,
        **kwargs
    ):
        """Plot a :py:class:`~radis.spectrum.spectrum.Spectrum` object.

        Parameters
        ----------
        var: variable (`absorbance`, `transmittance`, `transmittance_noslit`, etc.)
            For full list see :py:meth:`~radis.spectrum.spectrum.Spectrum.get_vars()`.
            If ``None``, plot the first thing in the Spectrum. Default ``None``.
        wunit: ``'default'``, ``'nm'``, ``'cm-1'``, ``'nm_vac'``,
            wavelength air, wavenumber, or wavelength vacuum. If ``'default'``,
            Spectrum :py:meth:`~radis.spectrum.spectrum.Spectrum.get_waveunit` is used.
        Iunit: unit for variable
            if `default`, default unit for quantity `var` is used.
            for radiance, one can use per wavelength (~ `W/m2/sr/nm`) or
            per wavenumber (~ `W/m2/sr/cm-1`) units


        Other Parameters
        ----------------
        show_points: boolean
            show calculated points. Default ``True``.
        nfig: int, None, or 'same'
            plot on a particular figure. 'same' plots on current figure. For
            instance::

                s1.plot()
                s2.plot(nfig='same')

        show_medium: bool, ``'vacuum_only'``
            if ``True`` and ``wunit`` are wavelengths, plot the propagation medium
            in the xaxis label (``[air]`` or ``[vacuum]``). If ``'vacuum_only'``,
            plot only if ``wunit=='nm_vac'``. Default ``'vacuum_only'``
            (prevents from inadvertently plotting spectra with different propagation
            medium on the same graph).
        yscale: 'linear', 'log'
            plot yscale
        normalize: boolean,  or tuple.
            option to normalize quantity to 1 (ex: for radiance). Default ``False``
        plot_by_parts: bool
            if ``True``, look for discontinuities in the wavespace and plot
            the different parts without connecting lines. Useful for experimental
            spectra produced by overlapping step-and-glue. Additional parameters
            read from ``kwargs`` : ``split_threshold`` and ``cutwings``. See more in
            :py:func:`~radis.spectrum.utils.split_and_plot_by_parts`.
        force: bool
            plotting on an existing figure is forbidden if labels are not the
            same. Use ``force=True`` to ignore that.
        **kwargs: **dict
            kwargs forwarded as argument to plot (e.g: lineshape
            attributes: `lw=3, color='r'`)

        Returns
        -------
        line:
            line plot

        Examples
        --------
        Plot an :py:func:`~radis.spectrum.models.experimental_spectrum` in
        arbitrary units::

            s = experimental_spectrum(..., Iunit='mW/cm2/sr/nm')
            s.plot(Iunit='W/cm2/sr/cm-1')

        See more examples in :ref:`the plot Spectral quantities page <label_spectrum_plot>`.

        See Also
        --------
        :py:func:`~radis.spectrum.compare.plot_diff`,
        :ref:`the Spectrum page <label_spectrum>`
        """

        # Deprecated
        if "plot_medium" in kwargs:
            show_medium = kwargs.pop("plot_medium")
            warn(DeprecationWarning("`plot_medium` was renamed to `show_medium`"))

        # Check inputs, get defaults
        # ------

        if var in ["intensity", "intensity_noslit"]:
            raise ValueError("`intensity` not defined. Use `radiance` instead")

        if var is None:  # if nothing is defined, try these first:
            params = self.get_vars()
            if "radiance" in params:
                var = "radiance"
            elif "radiance_noslit" in params:
                var = "radiance_noslit"
            elif "transmittance" in params:
                var = "transmittance"
            elif "transmittance_noslit" in params:
                var = "transmittance_noslit"
            else:
                # or plot the first variable we find
                var = list(params)[0]
                if var.replace("_noslit", "") in params:  # favour convolved quantities
                    var = var.replace("_noslit", "")

        if wunit == "default":
            wunit = self.get_waveunit()
        wunit = cast_waveunit(wunit)
        # Get variable
        x, y = self.get(var, wunit=wunit, Iunit=Iunit)

        # Get labels
        xlabel = format_xlabel(wunit, show_medium)
        if Iunit == "default":
            try:
                Iunit0 = self.units[var]
            except KeyError:  # unit not defined in dictionary
                Iunit0 = "a.u"
            Iunit = Iunit0

        # cosmetic changes
        ylabel = "{0} ({1})".format(make_up(var), make_up_unit(Iunit, var))
        # Plot
        # -------
        if normalize:
            if isinstance(normalize, tuple):
                from radis.misc.arrays import norm_on

                wmin, wmax = normalize
                y = norm_on(y, x, wmin=wmin, wmax=wmax)
            else:
                # y /= y.max()
                y /= np.nanmax(y)
            Iunit = "norm"

        set_style("origin")
        if nfig == "same":
            nfig = plt.gcf().number
        fig = plt.figure(nfig)

        # If figure exist, ensures xlabel and ylabel are the same (prevents some
        # users errors if plotting difference units!)... Note that since
        # 'radiance' and 'radiance_noslit' are now plotted under the same name,
        # they cannot be differenced. But at least this allows user to plot
        # both on the same figure if they want to compare [and have the same unit]

        def clean_error_msg(string):
            string = string.replace(r"$^\mathregular{", "^")
            string = string.replace(r"}$", "")
            return string

        if not force and (fig.gca().get_xlabel().lower() not in ["", xlabel.lower()]):
            raise ValueError(
                "Error while plotting {0}. Cannot plot ".format(var)
                + "on a same figure with different xlabel: {0}, {1}".format(
                    clean_error_msg(fig.gca().get_xlabel()), clean_error_msg(xlabel)
                )
                + "Use force=True if you really want to plot"
            )
        label1 = clean_error_msg(fig.gca().get_ylabel().lower())
        label2 = clean_error_msg(ylabel.lower())
        if not force and (label1 not in ["", label2]):
            raise ValueError(
                "Error while plotting {0}. Cannot plot ".format(var)
                + "on a same figure with different ylabel: \n{0}\n{1}".format(
                    clean_error_msg(fig.gca().get_ylabel()), clean_error_msg(ylabel)
                )
                + "\nUse force=True if you really want to plot"
            )

        # Add extra plotting parameters
        if "lw" not in kwargs and "linewidth" not in kwargs:
            kwargs["lw"] = 0.5
        # Add a label. Not shown by default but User can set it if using plt.legend()
        # (useful when plotting multiple plots on same figure)
        label = kwargs.pop("label", self.get_name())

        # Actual plot :
        # ... note: '-k' by default with style origin for first plot
        if not plot_by_parts:
            (line,) = plt.plot(x, y, label=label, **kwargs)
        else:
            (line,) = split_and_plot_by_parts(x, y, ax=fig.gca(), label=label, **kwargs)
            # note: split_and_plot_by_parts pops 'cutwing' & 'split_threshold' from kwargs

        if show_points:
            plt.plot(x, y, "o", color="lightgrey", **kwargs)

        # Labels
        plt.ticklabel_format(useOffset=False, axis="x")
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

        plt.yscale(yscale)

        if "label" in kwargs:
            plt.legend()
        fix_style(str("origin"))
        plt.show()
        return line

    def get_populations(self, molecule=None, isotope=None, electronic_state=None):
        """Return populations that are featured in the spectrum, either as
        upper or lower levels.

        Parameters
        ----------
        molecule: str, or None
            if None, only one molecule must be defined. Else, an error is raised
        isotope: int, or None
            isotope number. if None, only one isotope must be defined. Else,
            an error is raised
        electronic_state: str
            if None, only one electronic state must be defined. Else, an error
            is raised


        Returns
        -------
        pandas dataframe of levels, where levels are the index,
        and 'Evib' and 'nvib' are featured

        Notes
        -----
        Structure:

            {molecule: {isotope: {electronic_state: {'vib': pandas Dataframe,    # (copy of) vib levels
                                                     'rovib': pandas Dataframe,  # (copy of) rovib levels
                                                     'Ia': float    # isotopic abundance
                                                     }}}}

        (If Spectrum generated with RADIS, structure should match that of
        SpectrumFactory.get_populations())
        """

        # Check inputs, get default values
        populations = self.populations
        if populations is None or populations == {}:
            raise ValueError("Populations not defined")
        if type(populations) != dict:
            raise TypeError("Method defined for populations as dictionary")
        if molecule is None:
            if len(list(populations.keys())) != 1:
                raise ValueError(
                    "Please choose which molecule among: {0}".format(
                        list(populations.keys())
                    )
                )
            molecule = list(populations.keys())[0]
        if isotope is None:
            if len(list(populations[molecule].keys())) != 1:
                raise ValueError(
                    "Please choose which isotope among: {0}".format(
                        list(populations[molecule].keys())
                    )
                )
            isotope = list(populations[molecule].keys())[0]
        if electronic_state is None:
            if len(list(populations[molecule][isotope].keys())) != 1:
                raise ValueError(
                    "Please choose which electronic state among: {0}".format(
                        list(populations[molecule][isotope].keys())
                    )
                )
            electronic_state = list(populations[molecule][isotope].keys())[0]

        if __debug__:
            printdbg(
                "get vib populations for {0}({1})[iso{2}]".format(
                    molecule, electronic_state, isotope
                )
            )

        # Return
        return populations[molecule][isotope][electronic_state]

    def get_vib_levels(
        self, molecule=None, isotope=None, electronic_state=None, first=None
    ):
        """Return vibrational levels in the spectrum (energies, populations)

        Parameters
        ----------
        molecule: str, or None
            if None, only one molecule must be defined. Else, an error is raised
        isotope: int, or None
            isotope number. if None, only one isotope must be defined. Else,
            an error is raised
        electronic_state: str
            if None, only one electronic state must be defined. Else, an error
            is raised
        first: int, or 'all' or None
            only show the first N levels. If None or 'all', all levels are shown


        Returns
        -------
        out: pandas DataFrame
                pandas dataframe of levels, where levels are the index,
                and 'Evib' and 'nvib' are featured
        """

        pops = self.get_populations(
            molecule=molecule, isotope=isotope, electronic_state=electronic_state
        )

        try:
            vib_pops = pops["vib"]
        except KeyError:
            raise KeyError(
                "Vibrational levels not defined in this Spectrum object. "
                + "If using RADIS, make sure you chose export_populations='vib'"
            )

        if first is not None:
            if not "nvib" in vib_pops:
                raise KeyError(
                    "Vibrational populations (nvib) not calculated in this "
                    + "Spectrum object. Cant get first most populated levels. "
                    + "If using RADIS, make sure you used a non_eq_spectrum "
                    + "calculation"
                )
            if first == "all":
                first = None
            out = vib_pops.sort_values(by="nvib", ascending=False)[:first]

        else:
            out = vib_pops

        # Return
        return out

    def get_rovib_levels(
        self, molecule=None, isotope=None, electronic_state=None, first=None
    ):
        """Return rovibrational levels calculated in the spectrum (energies,
        populations)

        Parameters
        ----------
        molecule: str, or None
            if None, only one molecule must be defined. Else, an error is raised
        isotope: int, or None
            isotope number. if None, only one isotope must be defined. Else,
            an error is raised
        electronic_state: str
            if None, only one electronic state must be defined. Else, an error
            is raised
        first: int, or 'all' or None
            only show the first N levels. If None or 'all', all levels are shown


        Returns
        -------
        out: pandas DataFrame
                pandas dataframe of levels, where levels are the index,
                and 'Evib' and 'nvib' are featured
        """

        pops = self.get_populations(
            molecule=molecule, isotope=isotope, electronic_state=electronic_state
        )

        try:
            rovib_pops = pops["rovib"]
        except KeyError:
            raise KeyError(
                "Vibrational levels not defined in this Spectrum object. "
                + "If using RADIS, make sure you chose export_populations='vib'"
            )

        if first is not None:
            if not "n" in rovib_pops:
                raise KeyError(
                    "Rovibrational populations (n) not calculated in this "
                    + "Spectrum object. Cant get first most populated levels. "
                    + "If using RADIS, make sure you used a non_eq_spectrum "
                    + "calculation"
                )
            if first == "all":
                first = None
            out = rovib_pops.sort_values(by="n", ascending=False)[:first]

        else:
            out = rovib_pops

        # Return
        return out

    def plot_populations(
        self, what=None, nunit="", correct_for_abundance=False, **kwargs
    ):
        """Plots vib populations if given and format is valid.

        Parameters
        ----------
        what: 'vib', 'rovib', None
            if None plot everything
        nunit: '', 'cm-3'
            plot either in a fraction of vibrational levels, or a molecule
            number in in cm-3
        correct_for_abundance: boolean
            if ``True``, multiplies each population by the isotopic abundance
            (as it is done during the calculation of emission integral)
        kwargs: **dict
            are forwarded to the plot
        """

        # Check input, get defaults
        pops = self.populations
        if not isinstance(pops, dict):
            raise ValueError("Populations not defined in given Spectrum")

        # Plot function
        def _plot_elec_state(what, df, state_name, Ia, fig):
            if what == "vib":
                E, n, g = df["Evib"], df["nvib"], df["gvib"]
                ylabel = "n/g$_{vib}$"
                title = "Vibrational populations"
            elif what == "rovib":
                E, n, g = df["E"], df["n"], df["gj"]
                ylabel = "n/(2J+1)"
                title = "Rovibrational populations"

            if correct_for_abundance:
                n = n * Ia

            if nunit == "cm-3":
                from radis.phys.constants import k_b

                try:
                    P_mbar = self.conditions["pressure_mbar"]  # mbar
                    T = self.conditions["Tgas"]
                    mfrac = self.conditions["mole_fraction"]
                except KeyError:
                    raise KeyError(
                        "P_mbar (pressure), T (Tgas) and n (mole_fraction) "
                        + "are needed to calculate total number density in (cm-3)"
                    )
                N = P_mbar * 1e2 / k_b / T * mfrac * 1e-6
                n = n * N
                unitlabel = " [cm-3]"
            elif nunit == "":
                unitlabel = " [fraction]"
            else:
                raise ValueError("Unknown unit: {0}".format(nunit))

            # Plot
            if fig is None:
                fig = plt.figure()
                plt.xlabel("Energy (cm-1)")
                plt.ylabel(ylabel + unitlabel)
                plt.yscale("log")
                plt.title(title)
            ax = fig.gca()
            ax.plot(E, n / g, "o", label=state_name, **kwargs)

            return fig

        # Initialize figures, styles
        fig_vib = None
        fig_rovib = None
        set_style("origin")

        # Loop over all molecules, all isotopes, all electronic states
        # Note that the below works for both dict and pandas dataframe

        for molecule, isotopes in pops.items():
            for isotope, elec_states in isotopes.items():
                for elec_state, content in elec_states.items():
                    state_name = "{0}({1})(iso{2})".format(
                        molecule, elec_state, isotope
                    )

                    Ia = None
                    if correct_for_abundance:
                        if "Ia" in list(content.keys()):
                            Ia = content["Ia"]
                        else:
                            raise KeyError(
                                "Ia: isotopic abundance not defined in "
                                + "Spectrum populations"
                            )

                    for k in content.keys():

                        if k == "rovib" and (what == k or what is None):
                            df = self.get_rovib_levels(
                                molecule, isotope, elec_state, first="all"
                            )  # sort by n + check is defined)
                            fig_rovib = _plot_elec_state(
                                "rovib", df, state_name, Ia, fig_rovib
                            )

                        if k == "vib" and (what == k or what is None):
                            df = self.get_vib_levels(
                                molecule, isotope, elec_state, first="all"
                            )  # sort by n + check is defined)
                            fig_vib = _plot_elec_state(
                                "vib", df, state_name, Ia, fig_vib
                            )

        # Update
        for fig in [fig_vib, fig_rovib]:
            if fig is not None:
                ax = fig.gca()
                ax.legend()
                fix_style("origin", ax)

    # %% ------------------ Instrumental Slit Function ---------------------

    def apply_slit(
        self,
        slit_function,
        unit="nm",
        shape="triangular",
        center_wavespace=None,
        norm_by="area",
        mode="valid",
        plot_slit=False,
        store=True,
        slit_dispersion=None,
        slit_dispersion_warning_threshold=0.01,
        auto_recenter_crop=True,
        verbose=True,
        *args,
        **kwargs
    ):
        """Apply an instrumental slit function to all quantities in Spectrum.
        Slit function can be generated with usual shapes (see ``shape=``) or
        imported from an experimental slit function (path to a text file or
        numpy array of shape n*2). Convoluted spectra are cut on the edge
        compared to non-convoluted spectra, to remove side effects. See
        ``mode=`` to change this behaviour.

        Warning with units: read about ``'unit'`` and ``'return_unit'`` parameters.


        Parameters
        ----------
        slit_function: float or str or array
            If ``float``:
                generate slit function with FWHM of slit function (in nm or
                cm-1 depending on ``unit=``)
            If ``.txt``:
                import experimental slit function from .txt file: format must be 2-columns with
                wavelengths and intensity (doesn't have to be normalized)
            If ``array``:
                format must be 2-columns with wavelengths and intensity (doesn't have to be normalized)
        unit: ``'nm'`` or ``'cm-1'``
            unit of slit_function (FWHM, or imported file)
        shape: ``'triangular'``, ``'trapezoidal'``, ``'gaussian'``, or any of :data:`~radis.tools.slit.SLIT_SHAPES`
            which shape to use when generating a slit. Will call,
             respectively, :func:`~radis.tools.slit.triangular_slit`,
             :func:`~radis.tools.slit.trapezoidal_slit`,
             :func:`~radis.tools.slit.gaussian_slit`. Default 'triangular'
        center_wavespace: float, or ``None``
            center of slit when generated (in unit). Not used if slit is imported.
        norm_by: ``'area'``, ``'max'``
            normalisation type:
            - ``'area'`` normalizes the slit function to an area
              of 1. It conserves energy, and keeps the same units.
            - ``'max'`` normalizes the slit function to a maximum of 1.
              The convoluted spectrum units change (they are
              multiplied by the spectrum waveunit, e.g: a radiance
              non convoluted in mW/cm2/sr/nm on a wavelength (nm).
              range will yield a convoluted radiance in mW/cm2/sr.
              Note that the slit is set to 1 in the Spectrum wavespace
              (i.e: a Spectrum calculated in cm-1 will have a slit
              set to 1 in cm-1).
            Default ``'area'``
        mode: ``'valid'``, ``'same'``
           ``'same'`` returns output of same length as initial spectra,
            but boundary effects are still visible. ``'valid'`` returns
            output of length len(spectra) - len(slit) + 1, for
            which lines outside of the calculated range have
            no impact. Default ``'valid'``.

        Other Parameters
        ----------------
        auto_recenter_crop: bool
            if ``True``, recenter slit and crop zeros on the side when importing
            an experimental slit. Default ``True``.
            See :func:`~radis.tools.slit.recenter_slit`, :func:`~radis.tools.slit.crop_slit`
        plot_slit: boolean
            if ``True``, plot slit
        store: boolean
            if ``True``, store slit in the Spectrum object so it can be retrieved with
            :meth:`~radis.spectrum.spectrum.Spectrum.get_slit` and plot with
            :meth:`~radis.spectrum.spectrum.Spectrum.plot_slit`. Default ``True``
        slit_dispersion: func of (lambda, in ``'nm'``), or ``None``
            spectrometer reciprocal function : dλ/dx(λ)   (in ``nm``)
            If not ``None``, then the slit_dispersion function is used to correct the
            slit function for the whole range. Can be important if slit function
            was measured far from the measured spectrum  (e.g: a slit function
            measured at 632.8 nm will look broader at 350 nm because the spectrometer
            dispersion is higher at 350 nm. Therefore it should be corrected)
            Default ``None``

            .. warning::
                slit dispersion function is assumed to be given in ``nm``
                if your spectrum is stored in ``cm-1`` the wavenumbers are
                converted to wavelengths before being forwarded to the dispersion
                function

            See :func:`~radis.test.tools.test_slit.test_auto_correct_dispersion`
            for an example of the slit dispersion effect.

            A Python implementation of the slit dispersion:

            >>> def f(lbd):
            >>>    return  w/(2*f)*(tan(Φ)+sqrt((2*d/m/(w*1e-9)*cos(Φ))^2-1))

            Theoretical / References:

            >>> dλ/dx ~ d/mf    # at first order
            >>> dλ/dx = w/(2*f)*(tan(Φ)+sqrt((2*d/m/(w)*cos(Φ))^2-1))  # cf

            with:

            - Φ: spectrometer angle (°)
            - f: focal length (mm)
            - m: order of dispersion
            - d: grooves spacing (mm)   = 1/gr  with gr in (gr/mm)

            See Laux 1999 "Experimental study and modeling of infrared air plasma
            radiation" for more information

        slit_dispersion_warning_threshold: float
            if not ``None``, check that slit dispersion is about constant (< ``threshold`` change)
            on the calculated range. Default 0.01 (1%). See :func:`~radis.tools.slit.offset_dilate_slit_function`
        *args, **kwargs
            are forwarded to slit generation or import function
        verbose: boolean
            print stuff
        energy_threshold: float
             tolerance fraction when resampling. Default ``1e-3`` (0.1%)
             If areas before and after resampling differ by
             more than that an error is raised.


        Notes
        -----

        Units:

        the slit function is first converted to the wavespace (wavelength/wavenumber)
        that the Spectrum is stored in, and applied to the spectral quantities
        in their native wavespace.

        Implementation:

        :func:`~radis.tools.slit.convolve_with_slit` is applied to
        all quantities in :meth:`~radis.spectrum.spectrum.Spectrum.get_vars`
        that ends with _noslit. Generate a triangular instrumental slit function
        (or any other shape depending of shape=) with base
        ``slit_function_base`` (Uses the central wavelength of the spectrum
        for the slit function generation)

        We deal with several special cases (which makes the code
        a little heavy, but the method very versatile):

        - when slit unit and spectrum unit arent the same
        - when spectrum is not evenly spaced


        Examples
        --------
        ::

            s.apply_slit(1.2, 'nm')

        To manually apply the slit to a particular quantity use::

            wavenum, quantity = s['quantity']
            s['convolved_quantity'] = convolve_slit(wavenum, quantity,
                slit_function_base)

        See :func:`~radis.tools.slit.convolve_with_slit` for more details on Units and Normalization

        The slit is made considering the "center wavelength" which is
        the mean wavelength of the full spectrum you are applying it to.


        See Also
        --------
        :func:`~radis.tools.slit.get_slit_function`,
        :func:`~radis.tools.slit.convolve_with_slit`,
        :ref:`the Spectrum page <label_spectrum>`
        """
        # TODO: add warning if FWHM >= wstep(spectrum)/5

        from radis.tools.slit import (
            cast_waveunit,
            convolve_with_slit,
            get_slit_function,
            normalize_slit,
            offset_dilate_slit_function,
            remove_boundary,
        )

        # Check inputs
        # ---------
        if "slit_unit" in kwargs:
            unit = kwargs.pop("slit_unit")
            warn(DeprecationWarning("slit_unit was renamed unit"))

        unit = cast_waveunit(unit)

        varlist = [k for k in self.get_vars() if k.endswith("_noslit")]

        if len(varlist) == 0:
            raise AssertionError(
                "No variables to apply slit on. Variable names "
                + "to be convolved should end with _noslit"
            )

        # Forward relevant inputs to convolution instead of slit function generation
        kwargsconvolve = {}
        for kw in ["slit_dispersion", "verbose"]:
            if kw in kwargs:
                kwargsconvolve.update({kw: kwargs.pop(kw)})

        # For non evenyly distributed cases we take the minimum wstep among the
        # spectral range (and resample in convolve_with_slit)
        # Note: by construction all variables now have the same wavespace
        w = self._q["wavespace"]  # non convoluted wavespace
        wstep = abs(diff(w)).min()
        assert wstep > 0
        waveunit = self.get_waveunit()

        if __debug__:
            printdbg(
                "apply_slit: {0} in {1}, center `{2}`{1}, applied in waveunit {3}".format(
                    slit_function, unit, center_wavespace, waveunit
                )
            )

        if center_wavespace is None:
            # center_wavespace should be ~ unit
            center_wavespace = w[len(w) // 2]  # w ~ waveunit
            if waveunit == "cm-1" and unit == "nm":
                center_wavespace = cm2nm(center_wavespace)  # wavenum > wavelen
            elif waveunit == "nm" and unit == "cm-1":
                center_wavespace = nm2cm(center_wavespace)  # wavelen > wavenum

        # Get slit once and for all (and convert the slit unit
        # to the Spectrum `waveunit` if wavespaces are different)
        # -------
        wslit0, Islit0 = get_slit_function(
            slit_function,
            unit=unit,
            norm_by=norm_by,
            shape=shape,
            center_wavespace=center_wavespace,
            return_unit=waveunit,
            wstep=wstep,
            auto_recenter_crop=auto_recenter_crop,
            verbose=verbose,
            plot=plot_slit,
            *args,
            **kwargs
        )

        # Check if dispersion is too large
        # ----
        if waveunit == "nm":
            w_nm = w
            wslit0_nm = wslit0
        else:
            w_nm = cm2nm(w)
            wslit0_nm = cm2nm(wslit0)
        if slit_dispersion is not None:
            # add space (wings) on the side of each slice. This space is cut
            # after convolution, removing side effects. Too much space and we'll
            # overlap too much and loose performance. Not enough and we'll have
            # artifacts on the jonction. We use the slit width to be conservative.
            wings = len(wslit0)
            # correct if the slit is to be interpolated (because wstep is different):
            wings *= max(1, abs(int(np.diff(wslit0).mean() / wstep)))
            slice_windows, wings_min, wings_max = _cut_slices(
                w_nm, slit_dispersion, wings=wings
            )
        else:
            slice_windows = [np.ones_like(w, dtype=np.bool)]

        # Create dictionary to store convolved
        I_conv_slices = {}
        for qns in varlist:
            # Convolve and store the output in a new variable name (quantity name minus `_noslit`)
            # Create if requireds

            q = qns[:-7]  # new name  (minus '_noslit')
            w_conv_slices = []
            I_conv_slices[q] = []

        # Loop over all waverange slices (needed if slit changes over the spectral range)
        for islice, slice_window in enumerate(slice_windows):

            # Scale slit
            if slit_dispersion is not None:
                # apply spectrometer linear dispersion function.
                # dont forget it has to be added in nm and not cm-1
                wslit, Islit = offset_dilate_slit_function(
                    wslit0_nm,
                    Islit0,
                    w_nm[slice_window],
                    slit_dispersion,
                    threshold=slit_dispersion_warning_threshold,
                    verbose=verbose,
                )
                # Convert it back if needed
                if waveunit == "cm-1":
                    wslit = nm2cm(wslit)
                # We need to renormalize now that Islit has changed
                wslit, Islit = normalize_slit(wslit, Islit, norm_by=norm_by)
            else:
                wslit = wslit0
                Islit = Islit0  # no need to renormalize it

            # Apply to all variables
            # ---------
            for i, q in enumerate(I_conv_slices.keys()):
                # Convolve and store the output in a new variable name (quantity name minus `_noslit`)
                # Create if requireds

                qns = q + "_noslit"
                w_window = w[slice_window]
                I_window = self._q[qns][slice_window]

                # Apply convolution
                w_conv_window, I_conv_window = convolve_with_slit(
                    w_window,
                    I_window,
                    wslit,
                    Islit,
                    norm_by=None,  # already norm.
                    mode="same",  # dont loose information yet
                    waveunit=waveunit,
                    verbose=verbose,
                    assert_evenly_spaced=False,
                    # assumes Spectrum is correct by construction
                    **kwargsconvolve
                )

                # Crop wings to remove overlaps (before merging)
                if slit_dispersion is not None:
                    w_conv_window, I_conv_window = remove_boundary(
                        w_conv_window,
                        I_conv_window,
                        "crop",
                        crop_left=wings_min[islice],
                        crop_right=wings_max[islice],
                    )

                if i == 0:
                    w_conv_slices.append(w_conv_window)
                I_conv_slices[q].append(I_conv_window)

        # Merge and store all variables
        # ---------
        for q in I_conv_slices.keys():
            qns = q + "_noslit"
            I_not_conv = self._q[qns]

            # Merge all slices
            w_conv = np.hstack(w_conv_slices)
            I_conv = np.hstack(I_conv_slices[q])

            #            assert all(sorted(w_conv) == w_conv)

            # Crop to remove boundary effects (after merging)
            # this uses the mode='valid', 'same' attribute
            # ... if slit was imported, it has been interpolated on the Spectrum
            # ... grid and its initial length has changed: get the scaling factor
            # ... to remove the correct number of non valid points on the side
            scale_factor = abs((wslit0_nm[1] - wslit0_nm[0]) / (w_nm[1] - w_nm[0]))
            w_conv, I_conv = remove_boundary(
                w_conv,
                I_conv,
                mode,
                len_I=len(I_not_conv),
                len_I_slit_interp=int(len(Islit0) * scale_factor) + 1,
            )
            # Store
            self._q_conv["wavespace"] = w_conv
            self._q_conv[q] = I_conv

            # Get units
            if norm_by == "area":
                self.units[q] = self.units[qns]
            elif norm_by == "max":
                new_unit = (Unit(unit) * Unit(self.units[qns])).to_string()
                # because it's like if we multiplied by slit FWHM in the wavespace
                # it was generated
                self.units[q] = new_unit
            # Note: there was another mode called 'max2' where, unlike 'max',
            # unit was multiplied by [unit] not [return_unit]
            # Removed for simplification. You should stay with norm_by='area' anyway
            else:
                raise ValueError("Unknown normalization type: {0}".format(norm_by))

        # Sort if needed (sorting can be broken after applying corrected slits
        # on different slices)
        # | @EP: deactivated for the moment. It's dangerous to reorder because
        # | it creates features that could be confused with spectral features.
        #        if not is_sorted(w_conv) or not is_sorted_backward(w_conv):
        #            b = np.argsort(w_conv)
        #            for q in list(self._q_conv.keys()):
        #                self._q_conv[q] = self._q_conv[q][b]

        # Store slit in Spectrum, in the Spectrum unit
        if store:
            self._slit["wavespace"] = wslit0  # in 'waveunit'
            self._slit["intensity"] = Islit0

        # Update conditions
        self.conditions["slit_function"] = slit_function
        self.conditions["slit_unit"] = unit  # input slit unit
        self.conditions["slit_dispersion"] = slit_dispersion
        # TODO: probably removed after Spectrum is stored.
        self.conditions["norm_by"] = norm_by

        return self  # to be able to chain: s.apply_slit().plot()

    def get_slit(self, unit="same"):
        """Get slit function that was applied to the Spectrum.

        Returns
        -------

        wslit, Islit: array
            slit function with wslit in Spectrum ``waveunit``. See
            :meth:`~radis.spectrum.spectrum.Spectrum.get_waveunit`
        """

        if not unit in ["same", self.get_waveunit()]:
            raise NotImplementedError(
                "Unit must be Spectrum waveunit: {0}".format(self.get_waveunit())
            )

        # Make sure that slit is stored already
        try:
            wslit = self._slit["wavespace"]  # stored in Spectrum waveunit
            Islit = self._slit["intensity"]
        except KeyError:
            raise KeyError(
                "Slit function not found in Spectrum "
                + "conditions. Have you used Spectrum.apply_slit "
                + "with store=True?"
            )

        return wslit, Islit

    def plot_slit(self, wunit=None):
        """Plot slit function that was applied to the Spectrum.

        If dispersion was used (see :meth:`~radis.spectrum.spectrum.Spectrum.apply_slit`)
        the different slits are built again and plotted too (dotted).

        Parameters
        ----------

        wunit: ``'nm'``, ``'cm-1'``, or ``None``
            plot slit in wavelength or wavenumber. If ``None``, use the unit
            the slit in which the slit function was given. Default ``None``

        Returns
        -------

        fix, ax: matplotlib objects
            figure and ax

        See Also
        --------

        :ref:`the Spectrum page <label_spectrum>`
        """

        from radis.tools.slit import (
            normalize_slit,
            offset_dilate_slit_function,
            plot_slit,
        )

        # Check inputs
        assert wunit in ["nm", "cm-1", "nm_vac", None]
        # @dev: note: wunit in 'nm_vac' also implemented for consistency,
        # although not mentionned in docs.
        if wunit is None:
            wunit = self.conditions["slit_unit"]

        # Get slit arrays (in Spectrum.waveunit)
        wslit0, Islit0 = self.get_slit()  # as imported

        # Get slit unit
        norm_by = self.conditions["norm_by"]
        waveunit = self.get_waveunit()
        if norm_by == "area":
            Iunit = "1/{0}".format(waveunit)
        elif norm_by == "max":  # set maximum to 1
            Iunit = ""
        elif norm_by is None:
            Iunit = None
        else:
            raise ValueError(
                "Unknown normalization type: `norm_by` = {0}".format(norm_by)
            )

        # Plot in correct unit  (plot_slit deals with the conversion if needed)
        fig, ax = plot_slit(
            wslit0, Islit0, waveunit=waveunit, plot_unit=wunit, Iunit=Iunit
        )

        # Plot other slit functions if dispersion was applied:
        if "slit_dispersion" in self.conditions:
            slit_dispersion = self.conditions["slit_dispersion"]
            if slit_dispersion is not None:
                waveunit = self.get_waveunit()
                # Get slit in air wavelength:
                if waveunit == "nm":
                    wslit0_nm = wslit0
                elif waveunit == "nm_vac":
                    wslit0_nm = vacuum2air(wslit0)
                else:
                    wslit0_nm = cm2nm(wslit0)
                w_nm = self.get_wavelength(medium="air", which="non_convoluted")
                wings = len(
                    wslit0
                )  # note: hardcoded. Make sure it's the same as in apply_slit
                wings *= max(
                    1, abs(int(np.diff(wslit0_nm).mean() / np.diff(w_nm).mean()))
                )
                slice_windows, wings_min, wings_max = _cut_slices(
                    w_nm, slit_dispersion, wings=wings
                )

                # Loop over all waverange slices (needed if slit changes over the spectral range)
                for islice, slice_window in enumerate(slice_windows):

                    w_nm_sliced = w_nm[slice_window]
                    w_min = w_nm_sliced.min()
                    w_max = w_nm_sliced.max()

                    # apply spectrometer linear dispersion function.
                    # dont forget it has to be added in nm and not cm-1
                    wslit, Islit = offset_dilate_slit_function(
                        wslit0_nm,
                        Islit0,
                        w_nm[slice_window],
                        slit_dispersion,
                        threshold=0.01,
                        verbose=False,
                    )
                    # Convert it back if needed
                    if waveunit == "cm-1":
                        wslit = nm2cm(wslit)
                    elif waveunit == "nm_vac":
                        wslit = air2vacuum(wslit)
                    # We need to renormalize now that Islit has changed
                    wslit, Islit = normalize_slit(wslit, Islit, norm_by=norm_by)
                    plot_slit(
                        wslit,
                        Islit,
                        waveunit=waveunit,
                        plot_unit=wunit,
                        Iunit=Iunit,
                        ls="--",
                        title="Slit used on range {0:.2f}-{1:.2f} nm".format(
                            w_min, w_max
                        ),
                    )

        return fig, ax

    def line_survey(
        self,
        overlay=None,
        wunit="default",
        writefile=None,
        cutoff=None,
        *args,
        **kwargs
    ):
        """Plot Line Survey (all linestrengths used for calculation) Output in
        Plotly (html)

        Parameters
        ----------
        spec: Spectrum
            result from SpectrumFactory calculation (see spectrum.py)
        overlay: 'absorbance', 'transmittance', 'radiance', etc... or list of the above, or None
            overlay Linestrength with specified variable calculated in `spec`.
            Get the full list with the :meth:`~radis.spectrum.spectrum.Spectrum.get_vars`
            method. Default ``None``.
        wunit: ``'default'``, ``'nm'``, ``'cm-1'``, ``'nm_vac'``,
            wavelength air, wavenumber, or wavelength vacuum. If ``'default'``,
            Spectrum :py:meth:`~radis.spectrum.spectrum.Spectrum.get_waveunit` is used.
        medium: {'air', 'vacuum', 'default'}
            Choose whether wavelength are shown in air or vacuum. If ``'default'``
            lines are shown as stored in the spectrum.

        Other Parameters
        ----------------
        writefile: str
            if not ``None``, a valid filename to save the plot under .html format.
            If ``None``, use the ``fig`` object returned to show the plot.
        kwargs:: dict
            Other inputs are passed to :func:`~radis.tools.line_survey.LineSurvey`.
            Example below (see :py:func:`~radis.tools.line_survey.LineSurvey`
            documentation for more details):
        Iunit: `hitran`, `splot`
            Linestrength output units:

            - `hitran`: (cm-1/(molecule/cm-2))
            - `splot`: (cm-1/atm)   (Spectraplot units [2]_)

            Note: if not None, cutoff criteria is applied in this unit.
            Not used if plot is not 'S'

        barwidth: float
            With of bars in LineSurvey. Default 0.07



        Returns
        -------
        fig: a Plotly figure.
            If using a Jupyter Notebook, the plot will appear. Else, use ``writefile``
            to export to an html file.

        Plot in Plotly. See Output in [1]_


        Examples
        --------
        An example using the :class:`~radis.lbl.factory.SpectrumFactory` to generate a spectrum::

            from radis import SpectrumFactory
            sf = SpectrumFactory(
                                 wavenum_min=2380,
                                 wavenum_max=2400,
                                 mole_fraction=400e-6,
                                 path_length=100,  # cm
                                 isotope=[1],
                                 export_lines=True,    # required for LineSurvey!
                                 db_use_cached=True)
            sf.load_databank('HITRAN-CO2-TEST')
            s = sf.eq_spectrum(Tgas=1500)
            s.apply_slit(0.5)
            s.line_survey(overlay='radiance_noslit', barwidth=0.01)

        See the output in :ref:`Examples <label_examples>`


        References
        ----------
        .. [1] `RADIS Online Documentation (LineSurvey) <https://radis.readthedocs.io/en/latest/tools/line_survey.html>`__

        .. [2] `SpectraPlot <http://www.spectraplot.com/survey>`__


        See Also
        --------
        :func:`~radis.tools.line_survey.LineSurvey`,
        :ref:`the Spectrum page <label_spectrum>`
        """

        from radis.tools.line_survey import LineSurvey

        # Check inputs
        if wunit == "default":
            wunit = self.get_waveunit()

        def get_overlay(overlay):
            """Overlay line survey with a spectral quantity (like radiance or
            transmittance)"""

            if isinstance(overlay, str):  # either get it from the Spectrum
                if overlay not in self.get_vars():
                    raise AttributeError(
                        "{0} not in variables list: {1}".format(
                            overlay, self.get_vars()
                        )
                    )
                w, I = self.get(overlay, wunit=wunit)
                name = overlay
                units = self.units[overlay]
                return (w, I, name, units)

            else:  # Or use a given tuple or arrays
                try:
                    (w, I) = overlay
                except:
                    raise ValueError(
                        "Overlay has to be string, or (w,I) tuple of " + "arrays"
                    )
                return (w, I, "", "")

        if overlay is not None:

            if type(overlay) is not list:
                overlay = [overlay]

            overlay = [get_overlay(ov) for ov in overlay]

            return LineSurvey(
                self,
                overlay=overlay,
                wunit=wunit,
                writefile=writefile,
                cutoff=cutoff,
                *args,
                **kwargs
            )

        else:
            return LineSurvey(
                self, wunit=wunit, writefile=writefile, cutoff=cutoff, *args, **kwargs
            )

    def get_conditions(self):
        """Get all physical / computational parameters.

        See Also
        --------
        :py:meth:`~radis.spectrum.Spectrum.print_conditions`,
        :ref:`the Spectrum page <label_spectrum>`
        """

        return self.conditions

    def print_conditions(self, **kwargs):
        """Prints all physical / computational parameters. You can also simply
        print the Spectrum object directly::

            print(s)

        Parameters
        ----------
        kwargs: dict
            refer to :py:func:`~radis.spectrum.utils.print_conditions`

        See Also
        --------

        :py:meth:`~radis.spectrum.spectrum.Spectrum.get_conditions`,
        :py:func:`~radis.spectrum.utils.print_conditions`,
        :ref:`the Spectrum page <label_spectrum>`
        """

        return print_conditions(self.get_conditions(), self.cond_units, **kwargs)

    def store(
        self,
        path,
        discard=["lines", "populations"],
        compress=True,
        add_info=None,
        add_date=None,
        if_exists_then="error",
        verbose=True,
    ):
        """Save a Spectrum object in JSON format. Object can be recovered with
        :func:`~radis.tools.database.load_spec`. If many Spectrum are saved in a
        same folder you can view their properties with the :class:`~radis.tools.database.SpecDatabase`
        structure.

        Parameters
        ----------
        path: path to folder (database) or file
            if a folder, file is saved to database and name is generated automatically.
            if a file name, then Spectrum is saved to this file and the later
            formatting options dont apply
        file: str
            explicitely give a filename to save
        compress: boolean
            if ``False``, save under text format, readable with any editor.
            if ``True``, saves under binary format. Faster and takes less space.
            If ``2``, removes all quantities that can be regenerated with s.update(),
            e.g, transmittance if abscoeff and path length are given, radiance if
            emisscoeff and abscoeff are given in non-optically thin case, etc.
            Default ``True``.
        add_info: list
            append these parameters and their values if they are in conditions
            example::

                add_info = ['Tvib', 'Trot']

        discard: list of str
            parameters to exclude, for instance to save some memory for instance
            Default [`lines`, `populations`]: retrieved Spectrum looses the
            :meth:`~radis.spectrum.spectrum.Spectrum.line_survey` ability,
            and :meth:`~radis.spectrum.spectrum.Spectrum.plot_populations`
            (but it saves tons of memory!)
        if_exists_then: 'increment', 'replace', 'error'
            what to do if file already exists. If increment an incremental digit
            is added. If replace file is replaced (yeah). If error (or anything else)
            an error is raised. Default `error`


        Returns
        -------
        Returns filename used


        Notes
        -----
        If many spectra are stored in a folder, it may be time to set up a
        :class:`~radis.tools.database.SpecDatabase` structure to easily see all
        Spectrum conditions and get Spectrum that suits specific parameters.


        Implementation:

            Shouldnt rely on a Database. One may just want to store/load a Spectrum
            once.

        Examples
        --------
        Store a spectrum in compressed mode, regenerate quantities after loading::

            from radis import load_spec
            s.store('test.spec', compress=True)   # s is a Spectrum
            s2 = load_spec('test.spec')
            s2.update()                           # regenerate missing quantities


        See Also
        --------
        :class:`~radis.tools.database.SpecDatabase`,
        :func:`~radis.tools.database.load_spec`,
        :meth:`~radis.spectrum.spectrum.Spectrum.save`,
        :meth:`~radis.spectrum.spectrum.Spectrum.savetxt`

        """
        # TODO
        #
        # - in case a spectrometer linear dispersion function is used in
        # :meth:`~radis.spectrum.spectrum.Spectrum.apply_slit`, it probably isn't
        # stored with the current code. Find a workaround?

        from radis.tools.database import save

        if isinstance(discard, str):
            discard = [discard]

        return save(
            self,
            path,
            discard=discard,
            compress=compress,
            add_info=add_info,
            add_date=add_date,
            if_exists_then=if_exists_then,
            verbose=verbose,
        )

    def save(self, *args, **kwargs):
        """Alias to Spectrum.store.

        See Spectrum.store for documentation
        """

        return self.store(*args, **kwargs)

    def resample(
        self,
        w_new,
        unit="same",
        out_of_bounds="nan",
        if_conflict_drop="error",
        energy_threshold=1e-3,
        print_conservation=False,
        inplace=True,
        **kwargs
    ):
        """Resample spectrum over a new wavelength. Fills with transparent
        medium when out of bound (transmittance 1, radiance 0)

        .. warning::
            This may result in information loss. Resampling is done with
            oversampling and spline interpolation. These parameters can be adjusted,
            and energy conservation ensured with the appropriate parameters.

        Uses the :func:`radis.misc.signal.resample` function.


        Parameters
        ----------
        w_new: array,  or Spectrum
            new wavespace to resample the spectrum on. Must be inclosed in the
            current wavespace (we won't extrapolate)
            One can also give a Spectrum directly::

                s1.resample(s2.get_wavenumber())
                s1.resample(s2)            # also valid

        unit: ``'same'``, ``'nm'``, ``'cm-1'``, ``'nm_vac'``
            unit of new wavespace. It ``'same'`` it is assumed to be the current
            waveunit. Default ``'same'``. The spectrum waveunit is changed to this
            unit after resampling (i.e: a spectrum calculated and stored in `cm-1`
            but resampled in `nm` will be stored in `nm` from now on).
            If ``'nm'``, wavelength in air. If ``'nm_vac'``, wavelength in vacuum.
        out_of_bounds: 'transparent', 'nan', 'error'
            what to do if resampling is out of bounds. 'transparent': fills with
            transparent medium. 'nan': fill with nan. 'error': raises an error.
            Default 'nan'
        if_conflict_drop: 'error', 'convoluted', 'non_convoluted'
            There is a problem if both convoluted and non convoluted (*no_slit)
            quantities coexists, as they aren't scaled on the same wavespace
            grid. If 'error' an error is raised. If 'convoluted', convoluted
            quantities will be dropped. If 'non_convoluted' non convoluted quantities
            are dropped. Default 'error'
        medium: 'air', 'vacuum', or 'default'
            in which medium is the new waverange is calculated if it is given
            in 'nm'. Ignored if unit='cm-1'


        Other Parameters
        ----------------
        *Inputs forwarded to :func:`radis.misc.signal.resample`*

        energy_threshold: float
            if energy conservation (integrals) is above this threshold, raise an
            error.
        print_conservation: boolean
            if ``True``, prints energy conservation. Default ``False``.
        inplace: boolean
            if ``True``, modifies the Spectrum object directly. Else, returns
            a copy. Default ``True``.
        **kwargs: **dict
            all other arguments are sent to :func:`~radis.misc.signal.resample`

        Returns
        -------
        s: Spectrum
            resampled Spectrum object. If using ``inplace=True``, the Spectrum
            object has been modified anyway.

        See Also
        --------
        :func:`~radis.misc.signal.resample`
        """
        # TODO (but dangerous): reapply_slit at the end of the process if slit
        # is in conditions?

        if inplace:
            s = self
        else:
            s = self.copy()

        # Check inputs (check for deprecated)

        # ... see if convoluted / non convoluted values co-exist
        if "wavespace" in s._q and "wavespace" in s._q_conv:
            try:
                assert len(s._q["wavespace"]) == len(s._q_conv["wavespace"])
                assert np.allclose(s._q["wavespace"], s._q_conv["wavespace"])
            except AssertionError:  # wavespaces are not the same.
                if if_conflict_drop == "convoluted":
                    for q in list(s._q_conv.keys()):
                        del s._q_conv[q]
                elif if_conflict_drop == "non_convoluted":
                    for q in list(s._q.keys()):
                        del s._q[q]
                elif if_conflict_drop == "error":
                    raise ValueError(
                        "Cant resample as there are convoluted and non "
                        + "convoluted quantities in the Spectrum object (and "
                        + "wavespace are not the same). Use "
                        + "`if_conflict_drop='convoluted' or 'non_convoluted'`"
                    )
                else:
                    raise ValueError(
                        "Unknown value for if_conflict_drop: {0}".format(
                            if_conflict_drop
                        )
                    )

        # Get wavespace units
        stored_waveunit = s.get_waveunit()  # spectrum unit
        if unit == "same":  # resampled unit
            unit = stored_waveunit
        else:
            unit = cast_waveunit(unit)

        # Get output wavespace (it's the w_new array, unless w_new is a Spectrum)
        if isinstance(w_new, Spectrum):
            if unit == "nm":
                w_new = w_new.get_wavelength(medium="air")
            elif unit == "nm_vac":
                w_new = w_new.get_wavelength(medium="vacuum")
            elif unit == "cm-1":
                w_new = w_new.get_wavenumber()
            else:
                raise ValueError(unit)
        else:  # wavespace already given as array:
            w_new = w_new

        # Get current waverange in output unit   -> w
        if unit == "nm":
            w = s.get_wavelength(medium="air")
        elif unit == "nm_vac":
            w = s.get_wavelength(medium="vacuum")
        elif unit == "cm-1":
            w = s.get_wavenumber()
        else:
            raise ValueError("Unknown unit: {0}".format(unit))

        # Update stored_waveunit to new unit
        if unit != stored_waveunit:
            s.conditions["waveunit"] = unit

        # Get wavespace
        update_q = "wavespace" in s._q
        update_q_conv = "wavespace" in s._q_conv

        # Now let's resample
        def get_filling(variable):
            """Get out of bounds values for spectral quantity `variable`"""
            if out_of_bounds == "transparent":
                # Fill with optically transparent medium
                if variable in ["transmittance", "transmittance_noslit"]:
                    fill_with = 1
                else:
                    fill_with = 0
            elif out_of_bounds == "nan":
                fill_with = "nan"
            elif out_of_bounds == "error":
                fill_with = "error"
            else:
                raise ValueError(
                    "Unexpected value for out_of_bound: {0}".format(out_of_bounds)
                )
            return fill_with

        # There are different cases depending on the unit of w_new
        # ... Note @devs: we're looping over dictionaries directly rather than
        # ... using the (safer) .get() function because it's much faster (the
        # ... air2vacuum conversion in particular is quite slow, but has been
        # ... done once for all with get_wavelength() above )
        if update_q:
            for (k, I) in s._q.items():
                if k == "wavespace":
                    continue
                fill_with = get_filling(k)
                Inew = resample(
                    w,
                    I,
                    w_new,
                    ext=fill_with,
                    energy_threshold=energy_threshold,
                    print_conservation=False,
                    **kwargs
                )
                s._q[k] = Inew
            # update wavespace
            s._q["wavespace"] = w_new

        if update_q_conv:
            for (k, I) in s._q_conv.items():
                if k == "wavespace":
                    continue
                fill_with = get_filling(k)
                Inew = resample(
                    w,
                    I,
                    w_new,
                    ext=fill_with,
                    energy_threshold=energy_threshold,
                    print_conservation=False,
                    **kwargs
                )
                s._q_conv[k] = Inew
            # update wavespace
            s._q_conv["wavespace"] = w_new

        return s

    # %% ======================================================================
    # Semi public functions
    # ----------------
    # Access object properties (but don't manipulate the spectrum itself)
    # XXX =====================================================================

    def get_waveunit(self):
        """Returns whether this spectrum is defined in wavelength (nm) or
        wavenumber (cm-1)"""

        return self.conditions["waveunit"]

    def is_at_equilibrium(self, check="warn", verbose=False):
        """Returns whether this spectrum is at (thermal) equilibrium. Reads the
        ``thermal_equilibrium`` key in Spectrum conditions. It does not imply
        chemical equilibrium (mole fractions are still arbitrary)

        If they are defined, also check that the following assertions are True:

            Tvib = Trot = Tgas
            self_absorption = True
            overpopulation = None

        If they are not, still trust the value in Spectrum conditions, but raise
        a warning.

        Other Parameters
        ----------------
        check: ``'warn'``, ``'error'``, ``'ignore'``
            what to do if Spectrum conditions dont match the given equilibrium state:
            raise a warning, raise an error, or just ignore and dont even check.
            Default ``'warn'``.
        verbose: bool
            if ``True``, print why is the spectrum is not at equilibrium, if
            applicable.
        """

        conditions = self.conditions

        # Get value
        try:
            equilibrium = conditions["thermal_equilibrium"]
        except KeyError:
            raise KeyError(
                "We need to know if Spectrum is at equilibrium, but "
                + "`thermal_equilibrium` is not defined in conditions. Please add the "
                + "value manually with s.conditions['thermal_equilibrium']=..."
            )

        if check == "ignore":
            return equilibrium

        # Check output match the rest of the spectrum conditions
        try:
            assert conditions["Tgas"] != r"N/A"
            assert conditions["Tvib"] == conditions["Tgas"]
            assert conditions["Trot"] == conditions["Tgas"]
            if "overpopulation" in conditions:
                assert conditions["overpopulation"] is None
            assert conditions["self_absorption"]  # is True

            guess = True

        except AssertionError:
            guess = False
            if verbose:
                # Print which equilibrium test failed
                print("Spectrum not at equilibrium because the following test failed:")
                import sys
                import traceback

                _, _, tb = sys.exc_info()
                tb_info = traceback.extract_tb(sys.exc_info()[2])
                print(tb_info[-1][-1])
                # @dev: see https://stackoverflow.com/a/11587247/5622825
        except KeyError as err:
            warn(
                "Condition missing to know if spectrum is at equilibrium: {0}".format(
                    err
                )
            )
            guess = not equilibrium

        if equilibrium != guess:
            msg = (
                "Declared value of equilibrium ({0}) does not match the infered one ({1})".format(
                    equilibrium, guess
                )
                + ". Update your Spectrum conditions"
            )
            if check == "warn":
                warn(msg)
            elif check == "error":
                raise AssertionError(msg)

        return equilibrium

    def is_optically_thin(self):
        """Returns whether the spectrum is optically thin, based on the value
        on the self_absorption key in conditions.

        If not given, raises an error
        """

        try:
            return not self.conditions["self_absorption"]
        except KeyError:
            raise KeyError(
                "We need to know if Spectrum is optically thin, but "
                + "`self_absorption` is not defined in conditions. Please add the "
                + "value manually with s.conditions['self_absorption']=..."
            )

    def copy(self, copy_lines=True, quantity="all"):
        """Returns a copy of this Spectrum object (performs a smart deepcopy)

        Parameters
        ----------
        copy_lines: bool
            default ``True``
        quantity: 'all', or one of 'radiance_noslit', 'absorbance', etc.
            if not 'all', copy only one quantity. Default ``'all'``
        """
        try:
            return self.__copy__(copy_lines=copy_lines, quantity=quantity)
        except MemoryError:
            raise MemoryError(
                "during copy of Spectrum. If you don't need them, "
                + "droping lines before copying may save a lot of space: "
                + "del s.lines ; or, use copy_lines=False"
            )

    def __copy__(self, copy_lines=True, quantity="all"):
        """Generate a new spectrum object.

        Note: using deepcopy would work but then the Spectrum object would be pickled
        and unpickled again. It's a little faster here

        Parameters
        ----------
        copy_lines: bool
            default ``True``
        quantity: 'all', or one of 'radiance_noslit', 'absorbance', etc.
            if not 'all', copy only one quantity. Default ``'all'``

        Notes
        -----
        Performance:

            deepcopy: 3.32 ms
            initial : 10 ms
            no check test, optimised: 1.8 ms
            ... asymptote: without evenly spaced check, without copies: 1.84 ms
        """

        # quantities = {s:(v[0].copy(), v[1].copy()) for (s,v) in self.items()}  #╪ 1.8 ms
        #        quantities = dict(self.items())   # 912 ns, not a copy but no need as
        #                                        # Spectrum() recreates a copy anyway
        if quantity == "all":
            quantities = dict(self._get_items())
        else:
            #            assert quantity in CONVOLUTED_QUANTITIES+NON_CONVOLUTED_QUANTITIES
            #            if not quantity in self.get_vars():
            #                raise ValueError("Spectrum {0} has no quantity '{1}'. Got: {2}".format(
            #                        self.get_name(), quantity, self.get_vars()))
            quantities = {
                quantity: self.get(quantity, wunit=self.get_waveunit())
            }  # dict(self._get_items())

        try:
            units = deepcopy(self.units)
        except AttributeError:
            units = None

        conditions = deepcopy(self.conditions)  # 39 us
        try:
            cond_units = deepcopy(self.cond_units)
        except AttributeError:
            cond_units = None

        lines = None
        if copy_lines:
            try:
                # 143 us  (2 ms with deepcopy(lines))
                lines = self.lines.copy(deep=True)
            except AttributeError:
                pass

        try:
            populations = self.populations
        except AttributeError:
            populations = None

        waveunit = self.get_waveunit()  # 163 ns
        name = self.name

        # Generate copied Spectrum
        s = Spectrum(  # 1.51 ms
            quantities=quantities,
            conditions=conditions,
            cond_units=cond_units,
            populations=populations,
            lines=lines,
            units=units,
            waveunit=waveunit,
            name=name,
            warnings=False,  # saves about 3.5 ms on the Performance test object
        )

        # Add extra information

        # ... file name (if exists)
        s.file = self.file

        # ... slit information
        try:
            wslit, Islit = self.get_slit()
            s._slit["wavespace"] = wslit  # in 'waveunit'
            s._slit["intensity"] = Islit
        except KeyError:
            # no slit to apply
            pass

        return s

    def compare_with(
        self,
        other,
        spectra_only=False,
        plot=True,
        wunit="default",
        verbose=True,
        rtol=1e-5,
        ignore_nan=False,
        ignore_outliers=False,
        normalize=False,
        **kwargs
    ):
        """Compare Spectrum with another Spectrum object.

        Parameters
        ----------
        other: type Spectrum
            another Spectrum to compare with
        spectra_only: boolean, or str
            if ``True``, only compares spectral quantities (in the same waveunit)
            and not lines or conditions. If str, compare a particular quantity
            name. If False, compare everything (including lines and conditions
            and populations). Default ``False``
        plot: boolean
            if ``True``, use plot_diff to plot all quantities for the 2 spectra
            and the difference between them. Default ``True``.
        wunit: ``"nm"``, ``"cm-1"``, ``"default"``
            in which wavespace to compare (and plot). If ``"default"``, natural wavespace
            of first Spectrum is taken.
        rtol: float
            relative difference to use for spectral quantities comparison
        ignore_nan: boolean
            if ``True``, nans are ignored when comparing spectral quantities
        ignore_outliers: boolean, or float
            if not False, outliers are discarded. i.e, output is determined by::

                out = (~np.isclose(I, Ie, rtol=rtol, atol=0)).sum()/len(I) < ignore_outliers

        normalize: bool
            Normalize the spectra to be ploted

        Other Parameters
        ----------------
        kwargs: dict
            arguments are forwarded to :func:`~radis.spectrum.compare.plot_diff`

        Returns
        -------
        equals: boolean
            return True if spectra are equal (respective to tolerance defined by
            rtol and other input conditions)


        Examples
        --------
        Compare two Spectrum objects, or specifically the transmittance::

            s1.compare_with(s2)
            s1.compare_with(s2, 'transmittance')


        Note that you can also simply use `s1 == s2`, that uses
        :meth:`~radis.spectrum.spectrum.Spectrum.compare_with` internally::

            s1 == s2       # will return True or False


        See Also
        --------
        :func:`~radis.spectrum.compare.compare_spectra`
        """

        from radis.spectrum.compare import compare_spectra

        return compare_spectra(
            self,
            other,
            spectra_only=spectra_only,
            plot=plot,
            wunit=wunit,
            verbose=verbose,
            rtol=rtol,
            ignore_nan=ignore_nan,
            ignore_outliers=ignore_outliers,
            normalize=normalize,
            **kwargs
        )

    # %% ======================================================================
    # Private functions
    # ----------------
    # XXX =====================================================================

    def _init_annotations(self):
        """Annotations are used to give typing hints for get() and plot()
        functions, based on what quantities are available in the Spectrum
        object."""

        from radis.tools.slit import SLIT_SHAPES

        try:  # Python >3.6 only
            self.get.__annotations__["var"] = []
            self.plot.__annotations__["var"] = []
            self.apply_slit.__annotations__["shape"] = SLIT_SHAPES

        except AttributeError:
            pass  # old Python version

    def _add_quantity(self, name, w, I, warnings=True):
        """Add quantity.

        Note: creates a copy of the input array
        """

        assert len(w) == len(I)

        def check_wavespace(w):
            """If warnings, check that array is evenly spaced. Returns a copy
            of input array.

            Note: this check takes a lot of time!  (few ms)
            Is is not performed if warnings is False
            """
            if warnings:
                # Check Wavelength/wavenumber is evently spaced
                if not evenly_distributed(w, tolerance=1e-5):
                    warn(
                        "Wavespace is not evenly spaced ({0:.3f}%) for {1}.".format(
                            np.abs(np.diff(w)).max() / w.mean() * 100, name
                        )
                        + " This may create problems when convolving with slit function"
                    )
            return np.array(w)  # copy

        if name in CONVOLUTED_QUANTITIES:
            # Add wavespace
            if "wavespace" in self._q_conv:
                if warnings:
                    # Check new wavespace match the existing one
                    if not np.allclose(w, self._q_conv["wavespace"]):
                        raise ValueError(
                            "wavespace for {0} doesnt correspond to existing wavespace".format(
                                name
                            )
                            + " for convoluted quantities"
                        )
            else:
                self._q_conv["wavespace"] = np.array(w)  # copy
                # no need to check if wavespace is evenly spaced: we won't
                # apply the slit function again

            # Add quantity itself
            self._q_conv[name] = np.array(I)  # copy

        elif name in NON_CONVOLUTED_QUANTITIES:
            # Add wavespace
            if "wavespace" in self._q:
                if warnings:
                    # Check new wavespace match the existing one
                    if not np.allclose(w, self._q["wavespace"]):
                        raise ValueError(
                            "wavespace for {0} doesnt correspond to existing wavespace".format(
                                name
                            )
                            + " for non convoluted quantities"
                        )
            else:
                self._q["wavespace"] = check_wavespace(w)  # copy

            # Add quantity itself
            self._q[name] = np.array(I)  # copy

        else:
            raise ValueError(
                "Unknown quantity: {0}. Expected one of: {1}".format(
                    name, CONVOLUTED_QUANTITIES + NON_CONVOLUTED_QUANTITIES
                )
            )

        # also make the quantity accessible with s.[name] like Pandas dataframes (removed eventually)
        # setattr(self, name, quantity)   # Warning this makes another copy of it (it's a tuple!)

        # add to annotations   (Python >3.6)
        try:
            self.get.__annotations__["var"].append(name)
            self.plot.__annotations__["var"].append(name)
        except AttributeError:  # old Python version
            pass

    def __eq__(self, other):
        """Override the default Equals behavior."""
        return self.compare_with(other, verbose=False, plot=False)

    def __ne__(self, other):
        """Define a non-equality test."""
        return not self.__eq__(other)

    def __dir__(self):
        """Names shown with tab completion: remove certain attributes to
        simplify the use of this class (@minou)."""

        #        attrs = super(Spectrum, self).__dir__()
        attrs = dir(type(self))  # Python 2 and 3 compatible
        exclude = [
            "clear",
            "fromkeys",
            "items",
            "pop",
            "popitem",
            "setdefault",
            "values",
        ]

        return [k for k in attrs if not k in exclude]

    def __str__(self):
        """Print all Spectrum attributes."""

        # Print name
        print("Spectrum Name: ", self.get_name())

        # Print spectral quantities
        print("Spectral Quantities")
        print("-" * 40)
        for k, v in self._get_items().items():
            # print number of points with a comma separator
            print(
                " " * 2,
                k,
                "\t({0:,d} points{1})".format(
                    len(v[0]),
                    ", {0} nans".format(count_nans(v[1]))
                    if count_nans(v[1]) > 0
                    else "",
                ),
            )

        # Print populations
        print("Populations Stored")
        print("-" * 40)
        try:
            for k, v in self.populations.items():
                print(" " * 2, k, "\t\t", list(v.keys()))
        except:
            pass

        # Print conditions
        self.print_conditions()

        return ""  # self.print_conditions()

    def take(self, var):
        """
        Parameters
        ----------
        var : str
            spectral quantity

        Returns
        -------
        s: Spectrum
            same Spectrum with only the `var` spectral quantity

        Examples
        --------

        Use it to chain other commands ::

            s.take('radiance').normalize().plot()

        """

        return self.copy(quantity=var, copy_lines=True)

    # %% Add min, max, normalize operations

    def _get_unique_var(self, operation_name="algebraic"):
        quantities = self.get_vars()
        if len(quantities) > 1:
            raise KeyError(
                "There is an ambiguity with the Spectrum {0} operation. ".format(
                    operation_name
                )
                + "There should be only one var in Spectrum {0}. Got {1}\n".format(
                    self.get_name(), self.get_vars()
                )
                + "Use `s.take('transmittance')` or `s.take('radiance')`, etc. to extract the "
                "one spectral quantity you want."
            )
        elif len(quantities) == 0:
            raise KeyError(
                "No spectral quantity defined in Spectrum {0}".format(self.get_name())
            )
        else:
            var = quantities[0]
        return var

    def max(self):
        """Maximum of the Spectrum, if only one spectral quantity is
        available::

            s.max()

        Else, use :func:`~radis.spectrum.operations.Radiance`,
        :func:`~radis.spectrum.operations.Radiance_noslit`,
        :func:`~radis.spectrum.operations.Transmittance` or
        :func:`~radis.spectrum.operations.Transmittance_noslit`  ::

            Radiance(s).max()
        """

        var = self._get_unique_var(operation_name="max")
        w, I = self.get(var, wunit=self.get_waveunit(), copy=False)
        return I.max()

    def min(self):
        """Minimum of the Spectrum, if only one spectral quantity is available
        ::

            s.min()

        Else, use :func:`~radis.spectrum.operations.Radiance`,
        :func:`~radis.spectrum.operations.Radiance_noslit`,
        :func:`~radis.spectrum.operations.Transmittance` or
        :func:`~radis.spectrum.operations.Transmittance_noslit`  ::

            Radiance(s).min()
        """

        var = self._get_unique_var(operation_name="min")
        w, I = self.get(var, wunit=self.get_waveunit(), copy=False)
        return I.min()

    def normalize(
        self, normalize_how="max", wrange=(), wunit=None, inplace=False, force=False
    ):
        """Normalise the Spectrum, if only one spectral quantity is available.

        Parameters
        ----------
        normalize_how: ``'max'``, ``'area'``, ``'mean'``
            how to normalize. ``'max'`` is the default but may not be suited for very
            noisy experimental spectra. ``'area'`` will normalize the integral to 1.
            ``'mean'`` will normalize by the mean amplitude value
        wrange: tuple
            if not empty, normalize on this range
        wunit: ``"nm"``, ``"cm-1"``, ``"nm_vac"``
            unit of the normalisation range above. If ``None``, use the
            spectrum default waveunit.
        inplace: bool
            if ``True``, changes the Spectrum.

        Other Parameters
        ----------------
        force: boolean
            By default, normalizing some parametres such as transmittance
            is forbidden because considered non-physical. Use force=True
            if you really want to.

        Examples
        --------

            s.normalize("max", (4200, 4800), inplace=True)
        """

        from radis.spectrum.operations import multiply

        var = self._get_unique_var(operation_name="normalize")

        if var in ["transmittance", "transmittance_noslit"] and not force:
            raise ValueError(
                "Cannot normalize {0}. Use force=True if you really want.".format(var)
            )

        s = self

        if wunit is None:
            wunit = s.get_waveunit()

        if wrange is not None and len(wrange) > 0:
            wmin, wmax = wrange
            w, I = s.get(var, wunit=wunit, copy=False)  # (faster not to copy)
            b = (w > wmin) & (w < wmax)
            if normalize_how == "max":
                norm = np.nanmax(I[b])
                norm_unit = s.units[var]
            elif normalize_how == "mean":
                norm = np.nanmean(I[b])
                norm_unit = s.units[var]
            elif normalize_how == "area":
                norm = np.abs(nantrapz(I[b], w[b]))
                norm_unit = u.Unit(s.units[var]) * u.Unit(wunit)
            else:
                raise ValueError(
                    "Unexpected `normalize_how`: {0}".format(normalize_how)
                )

            out = multiply(s, 1 / norm, unit=norm_unit, inplace=inplace)

        else:
            if normalize_how == "max":
                norm = np.nanmax(s.get(var, copy=False)[1])
                norm_unit = s.units[var]

            elif normalize_how == "mean":
                norm = np.nanmean(s.get(var, copy=False)[1])
                norm_unit = s.units[var]

            elif normalize_how == "area":

                w, I = s.get(var, wunit=wunit, copy=False)
                norm = nantrapz(I, w)
                norm_unit = u.Unit(s.units[var]) * u.Unit(wunit)

            else:
                raise ValueError(
                    "Unexpected `normalize_how`: {0}".format(normalize_how)
                )
            # Ensure we use the same unit system!
            out = multiply(s, 1 / (norm * u.Unit(norm_unit)), inplace=inplace)

        return out

    # %% Define Spectrum Algebra
    # +, -, *, ^  operators

    # Notes:
    # s + s = 2*s    # valid only if optically thin

    # Other possibility:
    # use   s1>s2    for SerialSlabs
    # and   s1//s2   for MergeSlabs

    # Or:
    # use   s1*s2   for Serial Slabs
    # and   s1+s2   for MergeSlabs

    # For the moment the first version is implemented

    # Plus

    def __add__(self, other):
        """Override '+' behavior Add is defined as :

        - for numeric values: add a baseline (returns a copy)
        - for 2 Spectra: not defined (not physical)
        """

        if isinstance(other, float) or isinstance(other, int):
            from radis.spectrum.operations import add_constant

            return add_constant(self, other, inplace=False)
        elif isinstance(other, np.ndarray):
            from radis.spectrum.operations import add_array

            return add_array(self, other, inplace=True)
        elif isinstance(other, Spectrum):
            from radis.spectrum.operations import add_spectra

            return add_spectra(self, other)
        else:
            #            warn("You should'nt use the '+'. See '//' or '>' for more details", Warning)
            raise NotImplementedError(
                "+ not implemented for a Spectrum and a {0} object".format(type(other))
            )

    def __radd__(self, other):
        """Right side addition."""
        return self.__add__(other)

    def __iadd__(self, other):
        """Override '+=' behavior Add is defined as :

        - for numeric values: add a baseline (inplace)
        - for 2 Spectra: not defined (not physical)
        """
        if isinstance(other, float) or isinstance(other, int):
            from radis.spectrum.operations import add_constant

            return add_constant(self, other, inplace=True)
        elif isinstance(other, np.ndarray):
            from radis.spectrum.operations import add_array

            return add_array(self, other, inplace=True)
        else:
            warn("You should'nt use the '+'. See '//' or '>' for more details", Warning)
            raise NotImplementedError(
                "+ not implemented for a Spectrum and a {0} object".format(type(other))
            )

    # Minus

    def __sub__(self, other):
        """Override '-' behavior Add is defined as :

        - for numeric values: substract a baseline (returns a copy)
        - for 2 Spectra: defined only for baseline substraction
        """
        if isinstance(other, float) or isinstance(other, int):
            from radis.spectrum.operations import add_constant

            return add_constant(self, -other, inplace=False)
        elif isinstance(other, np.ndarray):
            from radis.spectrum.operations import add_array

            return add_array(self, -other, inplace=False)
        elif isinstance(other, Spectrum):
            from radis.spectrum.operations import substract_spectra

            return substract_spectra(self, other)
        else:
            raise NotImplementedError(
                "- not implemented for a Spectrum and a {0} object".format(type(other))
            )

    def __rsub__(self, other):
        """Right side substraction."""
        raise NotImplementedError(
            "right substraction (-) not implemented for Spectrum objects"
        )

    def __isub__(self, other):
        """Override '-=' behavior Add is defined as :

        - for numeric values: substract a baseline (inplace)
        - for 2 Spectra: defined only for baseline substraction
        """
        if isinstance(other, float) or isinstance(other, int):
            from radis.spectrum.operations import add_constant

            return add_constant(self, -other, inplace=True)
        elif isinstance(other, np.ndarray):
            from radis.spectrum.operations import add_array

            return add_array(self, -other, inplace=True)
        elif isinstance(other, Spectrum):
            from radis.spectrum.operations import substract_spectra

            return substract_spectra(self, other)
        else:
            raise NotImplementedError(
                "-= not implemented for a Spectrum and a {0} object".format(type(other))
            )

    # Times

    def __mul__(self, other):
        """Override '*' behavior Multiply is defined as :

        - for numeric values: multiply (equivalent to optically thin scaling)
          (only if in front, i.e:  2*s   works but s*2 is not implemented)
          (returns a copy)
        - for 2 Spectra: not defined
        """
        if (
            isinstance(other, float)
            or isinstance(other, int)
            or isinstance(other, u.quantity.Quantity)
        ):
            from radis.spectrum.operations import multiply

            return multiply(self, other, inplace=False)
        elif isinstance(other, Spectrum):
            raise NotImplementedError(
                "* not implemented for 2 Spectrum objects. Use > to combine them along the line of sight, as in SerialSlabs"
            )
        else:
            raise NotImplementedError(
                "* not implemented for a Spectrum and a {0} object".format(type(other))
            )

    def __rmul__(self, other):
        """Right side multiplication."""

        if (
            isinstance(other, float)
            or isinstance(other, int)
            or isinstance(other, u.quantity.Quantity)
        ):
            from radis.spectrum.operations import multiply

            return multiply(self, other, inplace=False)
        elif isinstance(other, Spectrum):
            raise NotImplementedError(
                "* not implemented for 2 Spectrum objects. Use > to combine them along the line of sight, as in SerialSlabs"
            )
        else:
            raise NotImplementedError(
                "right side * not implemented for a Spectrum and a {0} object".format(
                    type(other)
                )
            )

    def __imul__(self, other):
        """Override '*=' behavior Multiply is defined as :

        - for numeric values: multiply (equivalent to optically thin scaling)
          (only if in front, i.e:  s *= 2)  (modifies inplace)
        - for 2 Spectra: not defined
        """
        if (
            isinstance(other, float)
            or isinstance(other, int)
            or isinstance(other, u.quantity.Quantity)
        ):
            from radis.spectrum.operations import multiply

            return multiply(self, other, inplace=True)
        elif isinstance(other, Spectrum):
            raise NotImplementedError("* not implemented for 2 Spectrum objects. Use >")
        else:
            raise NotImplementedError(
                "*= not implemented for a Spectrum and a {0} object".format(type(other))
            )

    # Divide

    def __truediv__(self, other):
        """Override '/' behavior Divide is defined as :

        - for numeric values: divide algebrically (equivalent to optically thin scaling)
        """
        if isinstance(other, float) or isinstance(other, int):
            from radis.spectrum.operations import multiply

            return multiply(self, 1 / other, inplace=False)

        elif isinstance(other, u.quantity.Quantity):
            from radis.spectrum.operations import multiply

            return multiply(self, 1 / other.value, unit=1 / other.unit, inplace=False)

        else:
            raise NotImplementedError(
                "/ not implemented for a Spectrum and a {0} object".format(type(other))
            )

    def __rtruediv__(self, other):
        """Right side division."""

        raise NotImplementedError(
            "right side / not implemented for a Spectrum and a {0} object".format(
                type(other)
            )
        )

    def __itruediv__(self, other):
        """Override '/=' behavior Divide is defined as :

        - for numeric values: divide quantities algebrically
        (equivalent to optically thin scaling)
        """
        if isinstance(other, float) or isinstance(other, int):
            from radis.spectrum.operations import multiply

            return multiply(self, 1 / other, inplace=True)

        elif isinstance(other, u.quantity.Quantity):
            from radis.spectrum.operations import multiply

            return multiply(self, 1 / other.value, unit=1 / other.unit, inplace=True)

        else:
            raise NotImplementedError(
                "/= not implemented for a Spectrum and a {0} object".format(type(other))
            )

    # Line of sight operations

    def __bool__(self):
        """This prevents behaviors such as::

            s1 > s2 > s3

        which are actually interpreted by Python as "s1 > s2  and s2 > s3"
        and would return a wrong result  (probably s2 > s3  ? )
        """
        raise ArithmeticError(
            "A Spectrum cannot be evaluated as a boolean. "
            + "You may have tried using syntax such as `s1>s2>s3` "
            + "which Python interpets as `s1>s2 and s2>s3`. "
            + "Use `(s1>s2)>s3)` or SerialSlabs(s1, s2, s3) instead."
        )

    def __gt__(self, other):
        """Overloads '>' behavior no comparison: here we use > to define a
        ``Line of sight``.

        Examples
        --------
        s_plasma is seen through s_room::

            s = s_plasma > s_room

        Equivalent to::

            s = SerialSlabs(s_plasma, s_room)
        """
        if isinstance(other, Spectrum):
            from radis.los.slabs import SerialSlabs

            return SerialSlabs(self, other)
        else:
            raise NotImplementedError(
                "> not implemented for a Spectrum and a {0} object".format(type(other))
            )

    #    def __rgt__(self, other):
    #        ''' Right side > '''
    #
    #        if isinstance(other, Spectrum):
    #            from radis.los.slabs import SerialSlabs
    #            return SerialSlabs(other, self)
    #        else:
    #            raise NotImplementedError('right side > not implemented for a Spectrum and a {0} object'.format(
    #                    type(other)))

    def __floordiv__(self, other):
        """Overloads '//' behavior not a divison here: we use it to say that
        Slabs are ``in parallel``, i.e., as if their respective mole fractions
        were added in the same physical space.

        Won't work if they are not defined on the same waverange, but it's okay: let
        the user use MergeSlabs manually with the appropriate options

        Examples
        --------

        s_co2 added with s_co:

            s = s_co2 // s_co

        Equivalent to::

            s = MergeSlabs(s_co2, s_co)
        """

        if isinstance(other, Spectrum):
            from radis.los.slabs import MergeSlabs

            return MergeSlabs(self, other)
        else:
            raise NotImplementedError(
                "// not implemented for a Spectrum and a {0} object".format(type(other))
            )

    # the following is so that json_tricks.dumps and .loads can be used directly,
    # ie.:
    # >>> s == loads(s.dumps())
    # Does not work yet.
    # see https://github.com/mverleg/pyjson_tricks#class-instances
    #    def __json_encode__(self):
    #        ''' Called by json_tricks.dumps '''
    #        from radis.tools.database import _format_to_jsondict
    #        return _format_to_jsondict(self, discard=[], compress=[], verbose=True)
    #
    #    def __json_decode__(self, **attrs):
    #        from radis.tools.database import _json_to_spec
    #        print(attrs)
    #        raise
    #        return _json_to_spec(attrs)

    def __len__(self):
        """Length of a Spectrum object = length of the wavespace if unique,
        else raises an error"""

        # raises ValueError if both convolved and non convolved are defined
        try:
            return len(self._get_wavespace("any", copy=False))
        except ValueError:
            raise ValueError(
                "All quantities do not have the same length in the Spectrum : {0}".format(
                    {k: len(self.get(k)[0]) for k in self.get_vars()}
                )
            )


# %% Private functions

# to cut


def _cut_slices(w_nm, dispersion, threshold=0.01, wings=0):
    """used to cut a waverange into slices where dispersion does not very too
    much.

    Parameters
    ----------
    w_nm: numpy arrays
        wavelengths. If wavenumbers
    threshold: float
        must be a negative power of 10
    wings: int
        extend with that many points on each side of the slice. This space is cut
        after convolution (in apply_slit), removing side effects. Too much space and we'll
        overlap too much and loose performance. Not enough and we'll have
        artifacts on the jonction. A good number is to use the slit width, to be
        conservative.
    """

    # TODO: just test every 10 or 100
    blocs = dispersion(w_nm)
    diff = np.round(
        blocs[0] / blocs - 1, int(-np.log10(threshold))
    )  # difference in slit dispersion, +- 10%
    _, index_blocs = np.unique(diff, return_index=True)
    # check direction, add last element

    if index_blocs[0] < index_blocs[-1]:
        increment = 1
        index_blocs = np.hstack((index_blocs, len(w_nm)))
    else:
        increment = -1
        index_blocs = np.hstack((len(w_nm), index_blocs))
    #        index_blocs[-1] = None

    imins, imaxs = index_blocs[:-1].copy(), index_blocs[1:].copy()

    # Add wings
    if increment == 1:
        imins -= wings
        imaxs += wings
    else:
        imins += wings
        imaxs -= wings

    # add last if needed

    slices = []
    wings_min = []
    wings_max = []
    for imin, imax in zip(imins, imaxs):
        # Keep track of what was added on each side
        #        if imin <= - wings:
        #            wings_min.append(None)
        #            imin = None
        if imin <= 0:
            wings_min.append(wings + imin)
            imin = None
        else:
            wings_min.append(wings)
        #        if imax <= - wings:
        #            wings_max.append(None)
        #            imax = None
        if imax <= 0:
            wings_max.append(wings + imax)
            imax = None
        else:
            wings_max.append(wings)
        slice_w = np.zeros_like(w_nm, dtype=np.bool)
        slice_w[imin:imax:increment] = 1
        slices.append(slice_w)

    ##    # make sure we didnt miss anyone
    #    assert len(w_nm) == sum([slice_w.sum() for slice_w in slices])

    return slices[::increment], wings_min[::increment], wings_max[::increment]


# %% ======================================================================
# Test class function
# -------------------
# XXX =====================================================================

# Test class


def is_spectrum(a):
    """Returns whether a is a Spectrum object.

    Parameters
    ----------
    a: anything
        a Python object

    Returns
    -------
    True if a is a Spectrum object

    Notes
    -----
    is_spectrum compares the object class name (str): in some cases the Spectrum
    class gets imported twice (when databases are involved, mostly), and a purely
    isinstance() comparison fails
    """

    #    return isinstance(a, Spectrum)
    # removed: was used initially in the early RADIS development phase. Spectrum
    # object would not be recognized if the library was modified
    return isinstance(a, Spectrum) or repr(a.__class__) == repr(Spectrum)


# %% Test functions
if __name__ == "__main__":
    from radis.test.spectrum.test_spectrum import _run_testcases

    print("Test spectrum: ", _run_testcases(debug=False))
