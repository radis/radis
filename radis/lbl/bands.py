# -*- coding: utf-8 -*-
"""
Module that contains BandFactory, a class that features vibrational band-specific
functions, in particular the :meth:`~radis.lbl.bands.BandFactory.eq_bands` and
:meth:`~radis.lbl.bands.BandFactory.non_eq_bands` methods that return all
vibrational bands in a spectrum.

All of these are eventually integrated in SpectrumFactory which inherits from
BandFactory

Most methods are written in inherited class with the following inheritance scheme:

:py:class:`~radis.lbl.loader.DatabankLoader` > :py:class:`~radis.lbl.base.BaseFactory` >
:py:class:`~radis.lbl.broadening.BroadenFactory` > :py:class:`~radis.lbl.bands.BandFactory` >
:py:class:`~radis.lbl.factory.SpectrumFactory`

.. inheritance-diagram:: radis.lbl.factory.SpectrumFactory
   :parts: 1

Routine Listing
---------------

PUBLIC METHODS

- :meth:`~radis.lbl.bands.BandFactory.eq_bands`
- :meth:`~radis.lbl.bands.BandFactory.non_eq_bands`
- :meth:`~radis.lbl.bands.BandFactory.get_bands_weight`
- :meth:`~radis.lbl.bands.BandFactory.get_band_list`
- :meth:`~radis.lbl.bands.BandFactory.get_band`

PRIVATE METHODS

- _add_bands
- _broaden_lines_bands
- _broaden_lines_noneq_bands
- _calc_broadening_bands
- _calc_broadening_noneq_bands

----------


"""
# TODO: merge common parts of BandList.eq_bands  and SpectrumFactory.eq_spectrum,
# under a same function call

from time import time
from warnings import warn

import astropy.units as u
import numpy as np
import pandas as pd
from numpy import exp

from radis.api.hitranapi import HITRAN_CLASS1, get_molecule
from radis.lbl.broadening import BroadenFactory
from radis.lbl.labels import (
    vib_lvl_name_cdsd_pc,
    vib_lvl_name_cdsd_pcJN,
    vib_lvl_name_cdsd_pcN,
    vib_lvl_name_hitran_class1,
    vib_lvl_name_hitran_class5,
)

try:  # Proper import
    from .loader import KNOWN_DBFORMAT, KNOWN_LVLFORMAT
except ImportError:  # if ran from here
    from radis.lbl.loader import KNOWN_DBFORMAT, KNOWN_LVLFORMAT

from radis.misc.basics import all_in, is_float
from radis.misc.progress_bar import ProgressBar
from radis.misc.warning import reset_warnings
from radis.phys.constants import k_b
from radis.phys.units import convert_universal
from radis.phys.units_astropy import convert_and_strip_units
from radis.spectrum.equations import calc_radiance
from radis.spectrum.spectrum import Spectrum

# %% BandFactory


class BandFactory(BroadenFactory):
    """

    See Also
    --------

    :class:`~radis.lbl.factory.SpectrumFactory`
    """

    def __init__(self):

        super(BandFactory, self).__init__()

    # %% ======================================================================
    # PUBLIC METHODS
    # ------------------------
    # eq_bands            >>> get all calculation conditions
    # get_energy_levels   >>> return energy database
    #
    # =========================================================================

    # Vibrationaly specific calculations

    def eq_bands(
        self,
        Tgas,
        mole_fraction=None,
        diluent=None,
        path_length=None,
        pressure=None,
        levels="all",
        drop_lines=True,
    ):
        """Return all vibrational bands as a list of spectra for a spectrum
        calculated under equilibrium.

        By default, drop_lines is set to True so line_survey cannot be done on
        spectra. See drop_lines for more information

        Parameters
        ----------
        Tgas: float
            Gas temperature (K)
        mole_fraction: float
            database species mole fraction. If None, Factory mole fraction is used.
        path_length: float
            slab size (cm). If None, Factory mole fraction is used.
        pressure: float
            pressure (bar). If None, the default Factory pressure is used.

        Other Parameters
        ----------------
        levels: ``'all'``, int, list of str
            calculate only bands that feature certain levels. If ``'all'``, all
            bands are returned. If N (int), return bands for the first N levels
            (sorted by energy). If list of str, return for all levels in the list.
            The remaining levels are also calculated and returned merged together
            in the ``'others'`` key. Default ``'all'``
        drop_lines: boolean
            if False remove the line database from each bands. Helps save a lot
            of space, but line survey cannot be performed anymore. Default ``True``.

        Returns
        -------

        bands: list of :class:`~radis.spectrum.spectrum.Spectrum` objects

        Use .get(something) to get something among ['radiance', 'radiance_noslit',
                'absorbance', etc...]

        Or directly .plot(something) to plot it

        Notes
        -----

        Process:

        - Calculate line strenghts correcting the CDSD reference one.
        - Then call the main routine that sums over all lines
        """

        # Convert units
        Tgas = convert_and_strip_units(Tgas, u.K)
        path_length = convert_and_strip_units(path_length, u.cm)
        pressure = convert_and_strip_units(pressure, u.bar)

        # update defaults
        if path_length is not None:
            self.input.path_length = path_length
        if mole_fraction is not None:
            self.input.mole_fraction = mole_fraction
        if pressure is not None:
            self.input.pressure_mbar = pressure * 1e3
        if not is_float(Tgas):
            raise ValueError("Tgas should be a float or Astropy unit")
        assert type(levels) in [str, list, int]
        if type(levels) == str:
            assert levels == "all"
        # Temporary:
        if type(levels) == int:
            raise NotImplementedError

        # Get temperatures
        self.input.Tgas = Tgas
        self.input.Tvib = Tgas  # just for info
        self.input.Trot = Tgas  # just for info

        # Init variables
        pressure_mbar = self.input.pressure_mbar
        mole_fraction = self.input.mole_fraction
        path_length = self.input.path_length
        verbose = self.verbose

        # New Profiler object
        self._reset_profiler(verbose)

        # %% Retrieve from database if exists
        if self.autoretrievedatabase:
            s = self._retrieve_bands_from_database()
            if s is not None:
                return s

        # Print conditions
        if verbose:
            print("Calculating Equilibrium bands")
            self.print_conditions()

        # Start
        # %% Make sure database is loaded
        if self.df0 is None:
            raise AttributeError("Load databank first (.load_databank())")

        if not "band" in self.df0:
            self._add_bands()

        # %% Calculate the spectrum
        # ---------------------------------------------------

        self.profiler.start("band_calculation", 1)
        self._reinitialize()

        # --------------------------------------------------------------------

        # First calculate the linestrength at given temperature
        self.calc_linestrength_eq(Tgas)
        self._cutoff_linestrength()

        # ----------------------------------------------------------------------

        # Calculate line shift
        self.calc_lineshift()

        # ----------------------------------------------------------------------

        # Line broadening

        # ... generates molefraction for diluents
        self._generate_diluent_molefraction(mole_fraction, diluent)

        # ... calculate broadening  HWHM
        self._calc_broadening_HWHM()

        # ... generates all wstep related entities
        self._generate_wavenumber_arrays()

        # ... find weak lines and calculate semi-continuum (optional)
        I_continuum = self.calculate_pseudo_continuum()
        if I_continuum:
            raise NotImplementedError("pseudo continuum not implemented for bands")

        # ... apply lineshape and get absorption coefficient
        # ... (this is the performance bottleneck)
        wavenumber, abscoeff_v_bands = self._calc_broadening_bands()
        #    :         :
        #   cm-1    1/(#.cm-2)

        #            # ... add semi-continuum (optional)
        #            abscoeff_v_bands = self._add_pseudo_continuum(abscoeff_v_bands, I_continuum)

        # ----------------------------------------------------------------------
        # Remove certain bands
        if levels != "all":
            # Filter levels that feature the given energy levels. The rest
            # is stored in 'others'
            lines = self.df1
            # We need levels to be explicitely stated for given molecule
            assert hasattr(lines, "viblvl_u")
            assert hasattr(lines, "viblvl_l")
            # Get bands to remove
            merge_bands = []
            for (
                band
            ) in (
                abscoeff_v_bands
            ):  # note: could be vectorized with pandas str split. # TODO
                viblvl_l, viblvl_u = band.split("->")
                if not viblvl_l in levels and not viblvl_u in levels:
                    merge_bands.append(band)
            # Remove bands from bandlist and add them to `others`
            abscoeff_others = np.zeros_like(wavenumber)
            for band in merge_bands:
                abscoeff = abscoeff_v_bands.pop(band)
                abscoeff_others += abscoeff
            abscoeff_v_bands["others"] = abscoeff_others
            if verbose:
                print("{0} bands grouped under `others`".format(len(merge_bands)))

        # ----------------------------------------------------------------------
        # Generate spectra

        # Progress bar for spectra generation
        Nbands = len(abscoeff_v_bands)
        if self.verbose:
            print("Generating bands ({0})".format(Nbands))
        pb = ProgressBar(Nbands, active=self.verbose)
        if Nbands < 100:
            pb.set_active(False)  # hide for low line number

        # Generate spectra
        s_bands = {}
        for i, (band, abscoeff_v) in enumerate(abscoeff_v_bands.items()):

            # incorporate density of molecules (see equation (A.16) )
            density = mole_fraction * ((pressure_mbar * 100) / (k_b * Tgas)) * 1e-6
            #  :
            # (#/cm3)
            abscoeff = abscoeff_v * density  # cm-1

            # ==============================================================================
            # Warning
            # ---------
            # if the code is extended to multi-species, then density has to be added
            # before lineshape broadening (as it would not be constant for all species)
            # ==============================================================================

            # get absorbance (technically it's the optical depth `tau`,
            #                absorbance `A` being `A = tau/ln(10)` )
            absorbance = abscoeff * path_length

            # Generate output quantities
            transmittance_noslit = exp(-absorbance)
            emissivity_noslit = 1 - transmittance_noslit
            radiance_noslit = calc_radiance(
                wavenumber,
                emissivity_noslit,
                Tgas,
                unit=self.units["radiance_noslit"],
            )
            assert self.units["abscoeff"] == "cm-1"

            # ----------------------------- Export:

            lines = self.df1[self.df1.band == band]
            # if band == 'others': all lines will be None. # TODO
            populations = None  # self._get_vib_populations(lines)

            # Store results in Spectrum class
            if drop_lines:
                lines = None
                if self.save_memory:
                    try:
                        del self.df1  # saves some memory
                    except AttributeError:  # already deleted
                        pass
            conditions = self.get_conditions()
            # Add band name and hitran band name in conditions
            conditions.update({"band": band})

            if lines:

                def add_attr(attr):
                    if attr in lines:
                        if band == "others":
                            val = "N/A"
                        else:
                            # all have to be the same
                            val = lines[attr].iloc[0]
                        conditions.update({attr: val})

                add_attr("band_htrn")
                add_attr("viblvl_l")
                add_attr("viblvl_u")

            quantities = {
                "wavenumber": wavenumber,
                "abscoeff": abscoeff,
                "absorbance": absorbance,
                "emissivity_noslit": emissivity_noslit,
                "transmittance_noslit": transmittance_noslit,
                "radiance_noslit": radiance_noslit,
            }
            conditions["default_output_unit"] = self.input_wunit

            s = Spectrum(
                quantities=quantities,
                conditions=conditions,
                populations=populations,
                lines=lines,
                units=self.units,
                cond_units=self.cond_units,
                name=band,
                # dont check input (much faster, and Spectrum
                # is freshly baken so probably in a good format
                check_wavespace=False,
            )

            s_bands[band] = s

            pb.update(i)  # progress bar
        pb.done()

        self.profiler.stop("band_calculation", "Bands calculated")

        return s_bands

    def non_eq_bands(
        self,
        Tvib,
        Trot,
        Ttrans=None,
        mole_fraction=None,
        diluent=None,
        path_length=None,
        pressure=None,
        vib_distribution="boltzmann",
        rot_distribution="boltzmann",
        levels="all",
        return_lines=None,
    ):
        """Calculate vibrational bands in non-equilibrium case. Calculates
        absorption with broadened linestrength and emission with broadened
        Einstein coefficient.

        Parameters
        ----------
        Tvib: float
            vibrational temperature [K]
            can be a tuple of float for the special case of more-than-diatomic
            molecules (e.g: CO2)
        Trot: float
            rotational temperature [K]
        Ttrans: float
            translational temperature [K]. If None, translational temperature is
            taken as rotational temperature (valid at 1 atm for times above ~ 2ns
            which is the RT characteristic time)
        mole_fraction: float
            database species mole fraction. If None, Factory mole fraction is used.
        path_length: float
            slab size (cm). If None, Factory mole fraction is used.
        pressure: float
            pressure (bar). If None, the default Factory pressure is used.

        Other Parameters
        ----------------
        levels: ``'all'``, int, list of str
            calculate only bands that feature certain levels. If ``'all'``, all
            bands are returned. If N (int), return bands for the first N levels
            (sorted by energy). If list of str, return for all levels in the list.
            The remaining levels are also calculated and returned merged together
            in the ``'others'`` key. Default ``'all'``
        return_lines: boolean
            if ``True`` returns each band with its line database. Can produce big
            spectra! Default ``True``
            DEPRECATED. Now use export_lines attribute in Factory

        Returns
        -------

        Returns :class:`~radis.spectrum.spectrum.Spectrum` object

        Use .get(something) to get something among ['radiance', 'radiance_noslit',
        'absorbance', etc...]

        Or directly .plot(something) to plot it
        """

        # Convert units
        Tvib = convert_and_strip_units(Tvib, u.K)
        Trot = convert_and_strip_units(Trot, u.K)
        Ttrans = convert_and_strip_units(Ttrans, u.K)
        path_length = convert_and_strip_units(path_length, u.cm)
        pressure = convert_and_strip_units(pressure, u.bar)

        # check inputs, update defaults
        if path_length is not None:
            self.input.path_length = path_length
        if mole_fraction is not None:
            self.input.mole_fraction = mole_fraction
        if pressure is not None:
            self.input.pressure_mbar = pressure * 1e3
        if isinstance(Tvib, tuple):
            Tvib = tuple([convert_and_strip_units(T, u.K) for T in Tvib])
        elif not is_float(Tvib):
            raise TypeError(
                "Tvib should be float, or tuple (got {0})".format(type(Tvib))
            )
        singleTvibmode = is_float(Tvib)
        if not is_float(Trot):
            raise ValueError("Trot should be float.")
        assert type(levels) in [str, list, int]
        if type(levels) == str:
            assert levels == "all"
        else:
            if len(levels) != len(set(levels)):
                raise ValueError("levels list has duplicates")
        if not vib_distribution in ["boltzmann"]:
            raise ValueError("calculate per band not meaningful if not Boltzmann")
        # Temporary:
        if type(levels) == int:
            raise NotImplementedError
        if return_lines is not None:
            warn(
                DeprecationWarning(
                    "return_lines replaced with export_lines attribute in Factory"
                )
            )
            self.misc.export_lines = return_lines

        # Get translational temperature
        Tgas = Ttrans
        if Tgas is None:
            Tgas = Trot  # assuming Ttrans = Trot
        self.input.Tgas = Tgas
        self.input.Tvib = Tvib
        self.input.Trot = Trot

        # Init variables
        path_length = self.input.path_length
        mole_fraction = self.input.mole_fraction
        pressure_mbar = self.input.pressure_mbar
        verbose = self.verbose

        # New Profiler object
        self._reset_profiler(verbose)

        # %% Retrieve from database if exists
        if self.autoretrievedatabase:
            s = self._retrieve_bands_from_database()
            if s is not None:
                return s

        # Print conditions
        if verbose:
            print("Calculating Non-Equilibrium bands")
            self.print_conditions()

        self.profiler.start("band_calculation", 1)
        # %% Make sure database is loaded
        self._check_line_databank()
        self._calc_noneq_parameters(vib_distribution, singleTvibmode)

        if self.df0 is None:
            raise AttributeError("Load databank first (.load_databank())")

        if not "band" in self.df0:
            self._add_bands()

        # %% Calculate the spectrum
        # ---------------------------------------------------

        self._reinitialize()

        # ----------------------------------------------------------------------
        # Calculate Populations, Linestrength and Emission Integral
        # (Note: Emission Integral is non canonical quantity, equivalent to
        #  Linestrength for absorption)
        self.calc_populations_noneq(Tvib, Trot)
        self.calc_linestrength_noneq()
        self.calc_emission_integral()

        # ----------------------------------------------------------------------
        # Cutoff linestrength
        self._cutoff_linestrength()

        # ----------------------------------------------------------------------

        # Calculate lineshift
        self.calc_lineshift()

        # ----------------------------------------------------------------------

        # Line broadening

        # ... generates molefraction for diluents
        self._generate_diluent_molefraction(mole_fraction, diluent)

        # ... calculate broadening  HWHM
        self._calc_broadening_HWHM()

        # ... generates all wstep related entities
        self._generate_wavenumber_arrays()

        # ... find weak lines and calculate semi-continuum (optional)
        I_continuum = self.calculate_pseudo_continuum()
        if I_continuum:
            raise NotImplementedError("pseudo continuum not implemented for bands")

        # ... apply lineshape and get absorption coefficient
        # ... (this is the performance bottleneck)
        (
            wavenumber,
            abscoeff_v_bands,
            emisscoeff_v_bands,
        ) = self._calc_broadening_noneq_bands()
        #    :         :            :
        #   cm-1    1/(#.cm-2)   mW/sr/cm-1

        #            # ... add semi-continuum (optional)
        #            abscoeff_v_bands = self._add_pseudo_continuum(abscoeff_v_bands, I_continuum)

        # ----------------------------------------------------------------------
        # Remove bands
        if levels != "all":
            # Filter levels that feature the given energy levels. The rest
            # is stored in 'others'
            lines = self.df1
            # We need levels to be explicitely stated for given molecule
            assert hasattr(lines, "viblvl_u")
            assert hasattr(lines, "viblvl_l")
            # Get bands to remove
            merge_bands = []
            for (
                band
            ) in (
                abscoeff_v_bands
            ):  # note: could be vectorized with pandas str split. # TODO
                viblvl_l, viblvl_u = band.split("->")
                if not viblvl_l in levels and not viblvl_u in levels:
                    merge_bands.append(band)
            # Remove bands from bandlist and add them to `others`
            abscoeff_others = np.zeros_like(wavenumber)
            emisscoeff_others = np.zeros_like(wavenumber)
            for band in merge_bands:
                abscoeff = abscoeff_v_bands.pop(band)
                emisscoeff = emisscoeff_v_bands.pop(band)
                abscoeff_others += abscoeff
                emisscoeff_others += emisscoeff
            abscoeff_v_bands["others"] = abscoeff_others
            emisscoeff_v_bands["others"] = emisscoeff_others
            if verbose:
                print("{0} bands grouped under `others`".format(len(merge_bands)))

        # ----------------------------------------------------------------------
        # Generate spectra

        # Progress bar for spectra generation
        Nbands = len(abscoeff_v_bands)
        if self.verbose:
            print("Generating bands ({0})".format(Nbands))
        pb = ProgressBar(Nbands, active=self.verbose)
        if Nbands < 100:
            pb.set_active(False)  # hide for low line number

        # Create spectra
        s_bands = {}
        for i, band in enumerate(abscoeff_v_bands):
            abscoeff_v = abscoeff_v_bands[band]
            emisscoeff_v = emisscoeff_v_bands[band]

            # incorporate density of molecules (see equation (A.16) )
            density = mole_fraction * ((pressure_mbar * 100) / (k_b * Tgas)) * 1e-6
            #  :
            # (#/cm3)

            abscoeff = abscoeff_v * density  # cm-1
            emisscoeff = emisscoeff_v * density  # m/sr/cm3/cm-1

            # ==============================================================================
            # Warning
            # ---------
            # if the code is extended to multi-species, then density has to be added
            # before lineshape broadening (as it would not be constant for all species)
            # ==============================================================================

            # get absorbance (technically it's the optical depth `tau`,
            #                absorbance `A` being `A = tau/ln(10)` )

            # Generate output quantities
            absorbance = abscoeff * path_length  # (adim)
            transmittance_noslit = exp(-absorbance)

            if self.input.self_absorption:
                # Analytical output of computing RTE over a single slab of constant
                # emissivity and absorption coefficient
                b = abscoeff == 0  # optically thin mask
                radiance_noslit = np.zeros_like(emisscoeff)
                radiance_noslit[~b] = (
                    emisscoeff[~b] / abscoeff[~b] * (1 - transmittance_noslit[~b])
                )
                radiance_noslit[b] = emisscoeff[b] * path_length
            else:
                # Note that for k -> 0,
                radiance_noslit = emisscoeff * path_length  # (mW/cm2/sr/cm-1)

            # Convert `radiance_noslit` from (mW/sr/cm2/cm-1) to output unit
            radiance_noslit = convert_universal(
                radiance_noslit,
                from_unit="mW/cm2/sr/cm-1",
                to_unit=self.units["radiance_noslit"],
                wavenum=wavenumber,
                per_nm_is_like="mW/cm2/sr/nm",
                per_cm_is_like="mW/cm2/sr/cm-1",
            )
            # Convert 'emisscoeff' from (mW/sr/cm3/cm-1) to output unit
            emisscoeff = convert_universal(
                emisscoeff,
                from_unit="mW/cm3/sr/cm-1",
                to_unit=self.units["emisscoeff"],
                wavenum=wavenumber,
                per_nm_is_like="mW/cm3/sr/nm",
                per_cm_is_like="mW/cm3/sr/cm-1",
            )
            assert self.units["abscoeff"] == "cm-1"

            # ----------------------------- Export:

            lines = self.df1[self.df1.band == band]
            # Note: if band == 'others':  # for others: all will be None. # TODO. FIXME

            populations = self.get_populations(self.misc.export_populations)

            if not self.misc.export_lines:
                lines = None

            # Store results in Spectrum class
            if self.save_memory:
                try:
                    # saves some memory (note: only once 'lines' is discarded)
                    del self.df1
                except AttributeError:  # already deleted
                    pass
            conditions = self.get_conditions()
            conditions.update({"thermal_equilibrium": False})
            # Add band name and hitran band name in conditions

            def add_attr(attr):
                """# TODO: implement properly"""
                if lines is not None and attr in lines:
                    if band == "others":
                        val = "N/A"
                    else:
                        # all have to be the same
                        val = lines[attr].iloc[0]
                    conditions.update({attr: val})

            add_attr("band_htrn")
            add_attr("viblvl_l")
            add_attr("viblvl_u")

            quantities = {
                "wavenumber": wavenumber,
                "abscoeff": abscoeff,
                "absorbance": absorbance,
                "emisscoeff": emisscoeff,
                "transmittance_noslit": transmittance_noslit,
                "radiance_noslit": radiance_noslit,
            }
            conditions["default_output_unit"] = self.input_wunit

            s = Spectrum(
                quantities=quantities,
                conditions=conditions,
                populations=populations,
                lines=lines,
                units=self.units,
                cond_units=self.cond_units,
                name=band,
                # dont check input (much faster, and Spectrum
                # is freshly baken so probably in a good format
                check_wavespace=False,
            )

            s_bands[band] = s

            pb.update(i)  # progress bar
        pb.done()

        self.profiler.stop("band_calculation", "Bands calculated")

        return s_bands

    def get_bands_weight(self, showfirst=None, sortby="Ei"):
        """Show all bands by emission weight (fraction of total emission
        integral)

        Replace sortby with 'S' or 'int' to get different weights

        Note: this function works with .df1 (scaled lines) instead of .df0
        (reference lines) as we may want the weight for different temperatures.
        .df1 is updated after *eq_spectrum or non_eq_spectrum is called

        See Also
        --------

        :py:meth:`~radis.lbl.bands.BandFactory.get_bands`
        """

        # Check bands are labelled
        try:
            self.df1["band"]
        except KeyError:
            raise KeyError(
                "Please calculate bands first with SpectrumFactory.get_bands()"
            )

        df = self.df1
        dg = df.groupby("band")

        tot = df[sortby].sum()
        weight = (dg[sortby].sum() / tot).sort_values()[::-1]

        if showfirst is not None:
            return weight.head(showfirst)
        else:
            return weight

    def get_band_list(self):
        """Return all vibrational bands."""

        # Check bands are labelled
        try:
            self.df0["band"]
        except KeyError:
            self._add_bands()

        df = self.df0

        return df["band"].unique()

    def get_band(self, band):
        """Public function to get a particular vibrational band."""

        # Check bands are labelled
        try:
            self.df0["band"]
        except KeyError:
            self._add_bands()

        df = self.df0
        dg = df.groupby("band")

        return dg.get_group(band)

    # %% ======================================================================
    # PRIVATE METHODS
    # ---------------------------------
    # _add_bands
    # _broaden_lines_bands
    # _broaden_lines_noneq_bands
    # _calc_broadening_bands
    # _calc_broadening_noneq_bands

    # %% Line database functions: band specific

    def _add_bands(self):
        """Add a 'band' attribute for each line to allow parsing the lines by
        band with.

        >>> df0.groupby('band')


        Note on performance:
        -------

        test case (CDSD CO2 2380-2400 cm-1)
        - Initial: with .apply()   8.08 s ± 95.2 ms
        - with groupby(): 9s   worse!!
        - using simple (and more readable)    astype(str)  statements: 523 ms ± 19.6 ms
        """

        dbformat = self.params.dbformat
        lvlformat = self.params.levelsfmt
        verbose = self.verbose
        df = self.df0

        add_bands(df, dbformat, lvlformat, verbose=verbose)  # updates df

        return None  #  df already updated

    # Broadening functions: band specific

    def _broaden_lines_bands(self, df):
        """Divide over chuncks not to process to many lines in memory at the
        same time (note that this is not where the parallelisation is done: all
        lines are processed on the same core) Band specific version: returns a
        list of all broadened vibrational bands :

        Notes
        -----
        Implementation: there is no more splitting over line chuncks of given different
        (NotImplemented). This may result in large arrays and MemoryErrors for
        extreme spectral ranges. If that ever happens we may have to insert
        a chunck splitting loop in the band groupby loop

        See _calc_lineshape for more information
        """

        # Reactivate one-time warnings for new run
        reset_warnings(self.warnings)
        # --------------------------

        gb = df.groupby("band")

        abscoeff_bands = {}
        pb = ProgressBar(len(gb), active=self.verbose)
        optimization = self.params.optimization

        for i, (band, dg) in enumerate(gb):
            if optimization in ("simple", "min-RMS"):
                line_profile_LDM, wL, wG, wL_dat, wG_dat = self._calc_lineshape_LDM(dg)
                (wavenumber, absorption) = self._apply_lineshape_LDM(
                    dg.S.values,
                    line_profile_LDM,
                    dg.shiftwav.values,
                    wL,
                    wG,
                    wL_dat,
                    wG_dat,
                    optimization,
                )
            else:
                line_profile = self._calc_lineshape(dg)
                (wavenumber, absorption) = self._apply_lineshape(
                    dg.S.values, line_profile, dg.shiftwav.values
                )
            abscoeff_bands[band] = absorption
            pb.update(i)
        pb.done()

        return wavenumber, abscoeff_bands

    def _broaden_lines_noneq_bands(self, df):
        """Divide over chuncks not to process to many lines in memory at the
        same time (note that this is not where the parallelisation is done: all
        lines are processed on the same core) Band specific version: returns a
        list of all broadened vibrational bands.

        Notes
        -----
        Implementation: there is no more splitting over line chuncks of given different
        (NotImplemented). This may result in large arrays and MemoryErrors for
        extreme spectral ranges. If that ever happens we may have to insert
        a chunck splitting loop in the band groupby loop

        See _calc_lineshape for more information
        """

        # Reactivate one-time warnings for new run
        reset_warnings(self.warnings)
        # --------------------------

        abscoeff_bands = {}
        emisscoeff_bands = {}

        gb = df.groupby("band")
        optimization = self.params.optimization

        pb = ProgressBar(len(gb), active=self.verbose)
        for i, (band, dg) in enumerate(gb):
            if optimization in ("simple", "min-RMS"):
                line_profile_LDM, wL, wG, wL_dat, wG_dat = self._calc_lineshape_LDM(dg)
                (wavenumber, absorption) = self._apply_lineshape_LDM(
                    dg.S.values,
                    line_profile_LDM,
                    dg.shiftwav.values,
                    wL,
                    wG,
                    wL_dat,
                    wG_dat,
                    optimization,
                )
                (_, emission) = self._apply_lineshape_LDM(
                    dg.Ei.values,
                    line_profile_LDM,
                    dg.shiftwav.values,
                    wL,
                    wG,
                    wL_dat,
                    wG_dat,
                    optimization,
                )

            else:
                line_profile = self._calc_lineshape(dg)
                (wavenumber, absorption) = self._apply_lineshape(
                    dg.S.values, line_profile, dg.shiftwav.values
                )
                (_, emission) = self._apply_lineshape(
                    dg.Ei.values, line_profile, dg.shiftwav.values
                )
            abscoeff_bands[band] = absorption  #
            emisscoeff_bands[band] = emission
            pb.update(i)
        pb.done()

        return wavenumber, abscoeff_bands, emisscoeff_bands

    # %% Generate absorption profile which includes linebroadening factors

    def _calc_broadening_bands(self):
        """Loop over all lines, calculate lineshape, and returns the sum of
        absorption coefficient k=S*f over all lines. Band specific version:
        returns a list, one per vibrational band.

        For non-equilibrium, lineshape is calculated once and applied then
        to calculate absorption and emission coefficient.

        Returns
        -------

        abscoeff:  1/(#.cm-2)
            sum of all absorption coefficient k=1/(#.cm-2) for all lines in database
            `df` on the full calculation wavenumber range
        wavenumber: cm-1
            valid calculation wavenumber range

        Units
        -----

        `abscoeff` and `emisscoeff` still has to be multiplied by the total
        number density (cm-3) to get (cm-1/#) unit.

        """

        df = self.df1

        if self.verbose:
            print(("Now processing lines ({0})".format(len(df))))

        self.profiler.start("calc_broadening_eq_bands", 2)
        # Just some tests
        try:
            assert len(df.shape) == 2
        except AssertionError:
            warn(
                "Dataframe has only one line. Unexpected behaviour could occur"
                + " because Dataframes will be handled as Series and row/columns"
                + " may be inverted"
            )

        (wavenumber, abscoeff_bands) = self._broaden_lines_bands(df)

        self.profiler.stop("calc_broadening_eq_bands", "Calc Broadening Eq Bands")

        return wavenumber, abscoeff_bands

    def _calc_broadening_noneq_bands(self):
        """Loop over all lines, calculate lineshape, and returns the sum of
        absorption coefficient k=S*f over all lines. Band specific version:
        returns a list, one per vibrational band.

        For non-equilibrium, lineshape is calculated once and applied then
        to calculate absorption and emission coefficient.

        Returns
        -------

        wavenumber: cm-1
            full calculation wavenumber range
        abscoeff:  1/(#.cm-2)
            sum of all absorption coefficient k=1/(#.cm-2) for all lines in database
            `df` on the full calculation wavenumber range
        emisscoeff:  W/sr.cm
            sum of all broadened emission coefficients

        Units
        -----

        Both `abscoeff` and `emisscoeff` still have to be multiplied by the total
        number density (cm-3).
        """

        df = self.df1

        if self.verbose:
            print(("Now processing lines ({0})".format(len(df))))

        self.profiler.start("calc_broadening_noneq_bands", 2)
        # Just some tests
        try:
            assert len(df.shape) == 2
        except:
            warn(
                "Dataframe has only one line. Unexpected behaviour could occur"
                + " because Dataframes will be handled as Series and row/columns"
                + " may be inverted"
            )

        (
            wavenumber,
            abscoeff_bands,
            emisscoeff_bands,
        ) = self._broaden_lines_noneq_bands(df)

        self.profiler.stop(
            "calc_broadening_noneq_bands", "Calculate broadening noneq bands"
        )

        return wavenumber, abscoeff_bands, emisscoeff_bands

    def _retrieve_bands_from_database(self):
        """ """

        # TODO: options:
        # - Implement store / retrieve machinery for Bands
        # - The easiest way: generate a Database that contains a Band only. Assert
        #   all spectra have the same 'band' key and same conditions

        raise NotImplementedError("Retrieve bands from database not implemented")


# %% External functions


def docstring_parameter(*sub):
    def dec(obj):
        obj.__doc__ = obj.__doc__.format(*sub)
        return obj

    return dec


def add_bands(df, dbformat, lvlformat, verbose=True):
    """Assign all transitions to a vibrational band:

    Add 'band', 'viblvl_l' and 'viblvl_u' attributes for each line to allow
    parsing the lines by band with::

        df0.groupby('band')

    Parameters
    ----------

    df: pandas Dataframe
        Line (transitions) database
    dbformat: one of :data:`~radis.lbl.loader.KNOWN_DBFORMAT` : ``'cdsd```, ``'hitemp'``
        format of Line database
    lvlformat: 'cdsd`, 'hitemp'
        format of

    Returns
    -------

    None
        input Dataframe is updated inplace

    Examples
    --------

    Add transitions in a Dataframe based on CDSD (p, c, j, n) format::

        add_bands(df, 'cdsd-4000')

    Notes
    -----

    Performance with test case (CDSD CO2 2380-2400 cm-1):

    - Initial: with .apply()   8.08 s ± 95.2 ms
    - with groupby(): 9s   worse!!
    - using simple (and more readable)    astype(str)  statements: 523 ms ± 19.6 ms
    """

    # Check inputs
    if not dbformat in KNOWN_DBFORMAT:
        raise ValueError(
            "dbformat ({0}) should be one of: {1}".format(dbformat, KNOWN_DBFORMAT)
        )
    if not lvlformat in KNOWN_LVLFORMAT:
        raise ValueError(
            "lvlformat ({0}) should be one of: {1}".format(lvlformat, KNOWN_LVLFORMAT)
        )

    if verbose:
        t0 = time()
        print("... sorting lines by vibrational bands")

    # Calculate bands:
    if "id" in df:
        id = list(pd.unique(df["id"]))
        if len(id) > 1:
            raise ValueError(
                "Cant calculate vibrational bands for multiple " + "molecules yet"
            )  # although it's an easy fix. Just
            # groupby id
        molecule = get_molecule(id[0])

    else:
        id = df.attrs["id"]
        molecule = get_molecule(id)

    if molecule == "CO2":

        vib_lvl_name_hitran = vib_lvl_name_hitran_class5

        if lvlformat in ["cdsd-pc", "cdsd-pcN", "cdsd-hamil"]:

            # ensures that vib_lvl_name functions wont crash
            if dbformat not in [
                "cdsd-hitemp",
                "cdsd-4000",
                "hitran",
                "hitemp",
                "hitemp-radisdb",
            ]:
                raise NotImplementedError(
                    "lvlformat {0} not supported with dbformat {1}".format(
                        lvlformat, dbformat
                    )
                )

            # Use vibrational nomenclature of CDSD (p,c,j,n) or HITRAN (v1v2l2v3J)
            # depending on the Level Database.
            # In both cases, store the other one.

            # ... note: vib level in a CDSD (p,c,j,n) database is ambiguous.
            # ... a vibrational energy Evib can have been defined for every (p, c) group:
            if lvlformat in ["cdsd-pc"]:
                viblvl_l_cdsd = vib_lvl_name_cdsd_pc(df.polyl, df.wangl)
                viblvl_u_cdsd = vib_lvl_name_cdsd_pc(df.polyu, df.wangu)
            # ... or for every (p, c, N) group:
            elif lvlformat in ["cdsd-pcN"]:
                viblvl_l_cdsd = vib_lvl_name_cdsd_pcN(df.polyl, df.wangl, df.rankl)
                viblvl_u_cdsd = vib_lvl_name_cdsd_pcN(df.polyu, df.wangu, df.ranku)
            # ... or for every level (p, c, J ,N)  (that's the case if coupling terms
            # are used taken into account... it also takes a much longer time
            # to look up vibrational energies in the LineDatabase, warning!):
            elif lvlformat in ["cdsd-hamil"]:
                viblvl_l_cdsd = vib_lvl_name_cdsd_pcJN(
                    df.polyl, df.wangl, df.jl, df.rankl
                )
                viblvl_u_cdsd = vib_lvl_name_cdsd_pcJN(
                    df.polyu, df.wangu, df.ju, df.ranku
                )
            else:
                raise ValueError("Unexpected level format: {0}".format(lvlformat))

            band_cdsd = viblvl_l_cdsd + "->" + viblvl_u_cdsd

            df.loc[:, "viblvl_l"] = viblvl_l_cdsd
            df.loc[:, "viblvl_u"] = viblvl_u_cdsd
            df.loc[:, "band"] = band_cdsd

            # Calculate HITRAN format too (to store them))
            if all_in(["v1l", "v2l", "l2l", "v3l"], df):
                viblvl_l_hitran = vib_lvl_name_hitran(df.v1l, df.v2l, df.l2l, df.v3l)
                viblvl_u_hitran = vib_lvl_name_hitran(df.v1u, df.v2u, df.l2u, df.v3u)
                band_hitran = viblvl_l_hitran + "->" + viblvl_u_hitran

                df.loc[:, "viblvl_htrn_l"] = viblvl_l_hitran
                df.loc[:, "viblvl_htrn_u"] = viblvl_u_hitran
                df.loc[:, "band_htrn"] = band_hitran

        # 'radis' uses Dunham development based on v1v2l2v3 HITRAN convention
        elif lvlformat in ["radis"]:

            if dbformat not in [
                "hitran",
                "hitemp",
                "cdsd-hitemp",
                "cdsd-4000",
                "hitemp-radisdb",
            ]:
                raise NotImplementedError(
                    "lvlformat `{0}` not supported with dbformat `{1}`".format(
                        lvlformat, dbformat
                    )
                )

            # Calculate bands with HITRAN convention
            viblvl_l_hitran = vib_lvl_name_hitran(df.v1l, df.v2l, df.l2l, df.v3l)
            viblvl_u_hitran = vib_lvl_name_hitran(df.v1u, df.v2u, df.l2u, df.v3u)
            band_hitran = viblvl_l_hitran + "->" + viblvl_u_hitran

            df.loc[:, "viblvl_l"] = viblvl_l_hitran
            df.loc[:, "viblvl_u"] = viblvl_u_hitran
            df.loc[:, "band"] = band_hitran

        else:
            raise NotImplementedError(
                "Cant deal with lvlformat={0} for {1}".format(lvlformat, molecule)
            )

    elif molecule in HITRAN_CLASS1:  # includes 'CO'
        # Note. TODO. Move that in loader.py (or somewhere consistent with
        # classes defined in cdsd.py / hitran.py)

        if lvlformat in ["radis"]:

            # ensures that vib_lvl_name functions wont crash
            if dbformat not in ["hitran", "hitemp", "hitemp-radisdb"]:
                raise NotImplementedError(
                    "lvlformat {0} not supported with dbformat {1}".format(
                        lvlformat, dbformat
                    )
                )

            vib_lvl_name = vib_lvl_name_hitran_class1

            df.loc[:, "viblvl_l"] = vib_lvl_name(df["vl"])
            df.loc[:, "viblvl_u"] = vib_lvl_name(df["vu"])
            df.loc[:, "band"] = df["viblvl_l"] + "->" + df["viblvl_u"]

        else:
            raise NotImplementedError(
                "Lvlformat not defined for {0}: {1}".format(molecule, lvlformat)
            )

    else:
        raise NotImplementedError(
            "Vibrational bands not yet defined for molecule: "
            + "{0} with database format: {1}. ".format(molecule, dbformat)
            + "Update add_bands()"
        )

    if verbose:
        print(("... lines sorted in {0:.1f}s".format(time() - t0)))

    return


# %% Test


if __name__ == "__main__":

    from radis.test.lbl.test_bands import run_testcases

    print("test_bands.py:", run_testcases(plot=True))
