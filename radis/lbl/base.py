# -*- coding: utf-8 -*-
"""

Summary
-------

A class to aggregate methods to calculate spectroscopic parameter and
populations (and unload factory.py)

:py:class:`~radis.lbl.base.BaseFactory` is inherited by
:py:class:`~radis.lbl.broadening.BroadenFactory` eventually

Routine Listing
---------------


PUBLIC METHODS

- :py:meth:`radis.lbl.base.BaseFactory.print_conditions`         >>> get all calculation conditions
- :py:meth:`radis.lbl.base.BaseFactory.get_energy_levels`        >>> return energy database
- :py:meth:`radis.lbl.base.BaseFactory.plot_linestrength_hist`   >>>  plot distribution of linestrengths
- :py:meth:`radis.lbl.base.BaseFactory.plot_hist`                >>> same

PRIVATE METHODS - CALCULATE SPECTROSCOPIC PARAMETERS
(everything that doesnt depend on populations / temperatures)
(computation: work & update with 'df0' and called before eq_spectrum()  )

- :py:meth:`radis.lbl.base.BaseFactory._add_EvibErot`
- :py:meth:`radis.lbl.base.BaseFactory._add_EvibErot_CDSD`
- :py:meth:`radis.lbl.base.BaseFactory._add_EvibErot_RADIS_cls1`
- :py:meth:`radis.lbl.base.BaseFactory._add_Evib123Erot_RADIS_cls5`
- :py:meth:`radis.lbl.base.BaseFactory._add_ju`
- :py:meth:`radis.lbl.base.BaseFactory._add_Eu`
- :py:meth:`radis.lbl.base.BaseFactory._calc_noneq_parameters`
- :py:meth:`radis.lbl.base.BaseFactory.calc_weighted_trans_moment`
- :py:meth:`radis.lbl.base.BaseFactory.calc_einstein_coefficients`

PRIVATE METHODS - APPLY ENVIRONMENT PARAMETERS
(all functions that depends upon T or P)
(calculates populations, linestrength & radiance, lineshift)
(computation: work on df1, called by or after eq_spectrum() )

- :py:meth:`radis.lbl.base.BaseFactory.calc_lineshift`
- :py:meth:`radis.lbl.base.BaseFactory.calc_linestrength_eq`
- :py:meth:`radis.lbl.base.BaseFactory.calc_populations_eq`
- :py:meth:`radis.lbl.base.BaseFactory.calc_populations_noneq`
- :py:meth:`radis.lbl.base.BaseFactory.calc_linestrength_noneq`
- :py:meth:`radis.lbl.base.BaseFactory.calc_emission_integral`
- :py:meth:`radis.lbl.base.BaseFactory._cutoff_linestrength`

Most methods are written in inherited class with the following inheritance scheme:

:py:class:`~radis.lbl.loader.DatabankLoader` > :py:class:`~radis.lbl.base.BaseFactory` >
:py:class:`~radis.lbl.broadening.BroadenFactory` > :py:class:`~radis.lbl.bands.BandFactory` >
:py:class:`~radis.lbl.factory.SpectrumFactory`

.. inheritance-diagram:: radis.lbl.factory.SpectrumFactory
   :parts: 1


----------


"""
# TODO: move all CDSD dependant functions _add_Evib123Erot to a specific file for CO2.

import numpy as np
import pandas as pd
from astropy import units as u
from numpy import exp, pi
from psutil import virtual_memory

import radis

# TODO: rename in get_molecule_name
from radis.db.classes import get_molecule, get_molecule_identifier

try:  # Proper import
    from .loader import KNOWN_LVLFORMAT, DatabankLoader, df_metadata
except ImportError:  # if ran from here
    from radis.lbl.loader import KNOWN_LVLFORMAT, DatabankLoader, df_metadata

from radis.misc.arrays import anynan
from radis.misc.basics import all_in, is_float, transfer_metadata
from radis.misc.debug import printdbg
from radis.misc.log import printwarn
from radis.misc.plot import fix_style, set_style
from radis.misc.printer import printg
from radis.misc.utils import Default
from radis.misc.warning import OutOfBoundError
from radis.phys.constants import c_CGS, h_CGS, hc_k
from radis.phys.convert import cm2J, nm2cm, nm_air2cm
from radis.phys.units_astropy import convert_and_strip_units
from radis.spectrum.utils import print_conditions


class BaseFactory(DatabankLoader):

    # Output units
    units = {
        "absorbance": "",
        "abscoeff": "cm-1",
        "abscoeff_continuum": "cm-1",
        # TODO: deal with case where 'cm-1' is given as input for a Spectrum
        # (write a cast_unit of some kind)
        # different in Specair (mw/cm2/sr) because slit
        "radiance": "mW/cm2/sr/nm",
        # function is not normalized to conserve energy
        "radiance_noslit": "mW/cm2/sr/nm",  # it's actually a spectral radiance
        "emisscoeff": "mW/cm3/sr/nm",
        "emisscoeff_continuum": "mW/cm3/sr/nm",
        "emissivity": "",
        "emissivity_noslit": "",
        "transmittance": "",
        "transmittance_noslit": "",
    }

    # Calculation Conditions units
    cond_units = {
        "wavenum_min": "cm-1",
        "wavenum_max": "cm-1",
        "wavenum_min_calc": "cm-1",
        "wavenum_max_calc": "cm-1",
        "wstep": "cm-1",
        "wavelength_min": "nm",  # not defined as a variable. All calculations
        # are done with cm-1. Just for information
        "wavelength_max": "nm",
        "Tref": "K",
        "Tgas": "K",
        "Tvib": "K",
        "Trot": "K",
        "pressure_mbar": "mbar",
        "path_length": "cm",
        #        'slit_function_FWHM':   'nm',
        "cutoff": "cm-1/(#.cm-2)",
        "truncation": "cm-1",
        "neighbour_lines": "cm-1",
        # The later is never stored in Factory, but exported in Spectrum at the end of the calculation
        "calculation_time": "s",
    }

    def __init__(self):
        """

        .. inheritance-diagram:: radis.lbl.factory.SpectrumFactory
           :parts: 1

        """

        super(BaseFactory, self).__init__()  # initialize parent class

        # Define variable names
        # ... Note: defaults values are overwritten by SpectrumFactory input
        # ... values here are just to help autocompletion tools
        self.input.Tref = 296

        # all calculations are done in cm-1 in SpectrumFactory
        self.params.waveunit = "cm-1"
        # dont change this without making sure your line database
        # is correct, and the units conversions (ex: radiance)
        # are changed accordingly
        # Note that radiance are converted from ~ [mW/cm2/sr/cm-1]
        # to ~ [mW/cm/sr/nm]
        assert self.params.waveunit == self.cond_units["wstep"]

        self._export_continuum = False
        # private key to export abscoeff_continuum in the generated Spectrum.
        # so far continuum is not exported by default because rescale functions
        # are not defined yet. #TODO

    # %% ======================================================================
    # PUBLIC METHODS
    # ------------------------
    # print_conditions         >>> get all calculation conditions
    # get_energy_levels        >>> return energy database
    # plot_linestrength_hist   >>> get all linestrength distribution
    # plot_hist                >>> same
    #
    # =========================================================================

    def print_conditions(self, preprend=None):
        """Prints all physical / computational parameters. These are also
        stored in each result Spectrum.

        Parameters
        ----------

        preprend: str
            just to text to display before printing conditions
        """

        if preprend:
            print(preprend)

        conditions = self.get_conditions()

        return print_conditions(conditions, self.cond_units)

    def get_energy_levels(self, molecule, isotope, state, conditions=None):
        """Return energy levels database for given molecule > isotope > state
        (look up Factory.parsum_calc[molecule][iso][state])

        Parameters
        ----------

        molecule: str
            molecule name

        isotope: int
            isotope identifier

        state: str:
            electronic state

        conditions: str, or ``None``
            if not None, add conditions on which energies to retrieve, e.g:

            >>> 'j==0' or 'v1==0'

            Conditions are applied using Dataframe.query() method. In that case,
            ``get_energy_levels()`` returns a copy. Default ``None``

        Returns
        -------

        energies: pandas Dataframe
            a view of the energies stored in the Partition Function calculator
            for isotope iso. If conditions are applied, we get a copy

        See Also
        --------

        :meth:`~radis.lbl.base.BaseFactory.get_populations`
        """

        energies = self.get_partition_function_calculator(molecule, isotope, state).df

        if conditions is not None:
            energies = energies.query(conditions)

        return energies

    def plot_linestrength_hist(self, cutoff=None):
        """Plot linestrength distribution (to help determine a cutoff
        criteria)"""
        return self.plot_hist("df1", "S", axvline=np.log10(cutoff))

    def plot_hist(self, dataframe="df0", what="int", axvline=None):
        """Plot distribution of column ``what`` in ``dataframe``

        For instance, help determine a cutoff criteria ::

            plot_hist("df1", "int")

        Parameters
        ----------
        dataframe: 'df0', 'df1'
            which dataframe to plot (df0 is the loaded one, df1 the scaled one)
        what: str
            which feature to plot. Default ``'S'`` (scaled linestrength). Could also
            be ``'int'`` (reference linestrength intensity), ``'A'`` (Einstein coefficient)
        axvline: float
            if not ``None``, plot a vertical line at this position.
        """
        import matplotlib.pyplot as plt

        assert dataframe in ["df0", "df1"]
        plt.figure()
        df = getattr(self, dataframe)
        a = np.log10(np.array(df[what]))
        if np.isnan(a).any():
            printwarn("Nan values in log10(lines)")
        plt.hist(np.round(a[~np.isnan(a)]))
        if axvline is not None:
            plt.axvline(axvline, color="r")
        plt.xlabel("log10({0})".format(what))
        plt.ylabel("Count")
        plt.show()

    # %% ======================================================================
    # PRIVATE METHODS - CALCULATE SPECTROSCOPIC PARAMETERS
    # (everything that doesnt depend on populations / temperatures)
    # (computation: work & update with 'df0' and called before eq_spectrum()  )
    # ---------------------------------
    # _add_EvibErot
    # _add_EvibErot_CDSD
    # _add_EvibErot_RADIS_cls1
    # _add_Evib123Erot_RADIS_cls5
    # _add_ju
    # _add_Eu
    # _calc_noneq_parameters
    # calc_weighted_trans_moment
    # calc_einstein_coefficients
    # =========================================================================

    def assert_no_nan(self, df, column):
        """Assert there are no nan in the column.

        Crash with a nice explanation if one is found"""
        from radis.misc.printer import get_print_full

        try:
            assert not anynan(df[column])
        except AssertionError as err:
            index = np.isnan(df[column]).idxmax()
            if self.input.molecule == "CO2":
                fix_idea = (
                    "If using HITEMP2010 for CO2, some lines are unlabelled and therefore cannot be used at "
                    "equilibrium. This is a known issue of the HITEMP database and will soon be fixed in the "
                    "edition. In the meantime you can use:\n 'sf.df0.drop(sf.df0.index[sf.df0['v1u']==-1], inplace=True)' "
                    "where 'sf' is SpectrumFactory object"
                )
            raise AssertionError(
                "{0}=NaN in line database at index {1}".format(column, index)
                + " corresponding to Line:\n {0}".format(
                    get_print_full(df.loc[index]) + fix_idea
                )
            ) from err

    def _add_EvibErot(self, df, calc_Evib_harmonic_anharmonic=False):
        """Calculate Evib & Erot in Line dataframe.

        Parameters
        ----------

        df: DataFrame
            list of transitions

        Other Parameters
        ----------------

        calc_Evib_harmonic_anharmonic: boolean
            if ``True``, calculate harmonic and anharmonic components of
            vibrational energies (for Treanor distributions)
        """
        if self.verbose:
            print("Fetching Evib & Erot.")
            if self.verbose >= 2:
                printg(
                    "If using this code several"
                    + " times you should consider updating the database"
                    + " directly. See functions in factory.py "
                )

        from radis.db.classes import HITRAN_CLASS1, HITRAN_CLASS5

        # Different methods to get Evib and Erot:
        # fetch energies from precomputed CDSD levels: one Evib per (p, c) group
        if self.params.levelsfmt == "cdsd-pc":
            return self._add_EvibErot_CDSD_pc(
                df, calc_Evib_harmonic_anharmonic=calc_Evib_harmonic_anharmonic
            )

        # fetch energies from precomputed CDSD levels: one Evib per (p, c, N) group
        elif self.params.levelsfmt == "cdsd-pcN":
            return self._add_EvibErot_CDSD_pcN(
                df, calc_Evib_harmonic_anharmonic=calc_Evib_harmonic_anharmonic
            )

        # fetch energies from CDSD levels calculated from Hamiltonian: one Evib per (p, c, J, N) group
        # (that's necessary if Evib are different for all levels, which can be
        # the case if coupling terms are calculated)
        elif (
            self.params.levelsfmt == "cdsd-hamil"
        ):  # fetch energies from precomputed CDSD levels
            return self._add_EvibErot_CDSD_pcJN(
                df, calc_Evib_harmonic_anharmonic=calc_Evib_harmonic_anharmonic
            )

        # calculate directly with Dunham expansions, whose terms are included in
        # the radis.db database
        elif self.params.levelsfmt == "radis":
            molecule = self.input.molecule
            if molecule in HITRAN_CLASS1:  # class 1
                return self._add_EvibErot_RADIS_cls1(
                    df, calc_Evib_harmonic_anharmonic=calc_Evib_harmonic_anharmonic
                )
            elif molecule in HITRAN_CLASS5:
                if __debug__:
                    printdbg(
                        "placeholder: using getEvib123 while getEvib would be enough"
                    )
                # TODO: write simplified function: doesnt need to fetch Evib1,2,3 here
                # as only Evib is needed in this case
                if calc_Evib_harmonic_anharmonic:
                    return self._add_Evib123Erot_RADIS_cls5_harmonicanharmonic(df)
                else:
                    return self._add_Evib123Erot_RADIS_cls5(df)
            else:
                raise NotImplementedError(
                    "Molecules not implemented: {0}".format(molecule.name)
                )  # TODO

        else:
            raise NotImplementedError(
                "Impossible to calculate Evib Erot with given energy "
                + "format: {0}".format(self.params.levelsfmt)
            )

    def _add_Evib123Erot(self, df, calc_Evib_harmonic_anharmonic=False):
        """Calculate Evib & Erot in dataframe.

        Parameters
        ----------

        df: DataFrame
        """
        if self.verbose:
            print(
                "Fetching Evib & Erot. If using this code several"
                + " times you should consider updating the database"
                + " directly. See functions in factory.py "
            )

        from radis.db.classes import HITRAN_CLASS5

        if self.params.levelsfmt == "cdsd-pc":  # calculate from precomputed CDSD levels
            return self._add_Evib123Erot_CDSD_pc(df)

        elif (
            self.params.levelsfmt == "cdsd-pcN"
        ):  # calculate from precomputed CDSD levels
            raise NotImplementedError("3 Tvib mode for CDSD in pcN convention")  # TODO

        elif self.params.levelsfmt == "radis":  # calculate with Dunham expansions
            if self.input.molecule in HITRAN_CLASS5:  # class 5
                if calc_Evib_harmonic_anharmonic:
                    return self._add_Evib123Erot_RADIS_cls5_harmonicanharmonic(df)
                else:
                    return self._add_Evib123Erot_RADIS_cls5(df)
            else:
                raise NotImplementedError(
                    "Molecules other than HITRAN class 5 (CO2) not implemented"
                )  # TODO

        else:
            raise NotImplementedError(
                "Impossible to calculate Evib Erot with given energy "
                + "format: {0}".format(self.params.levelsfmt)
            )

    def _add_EvibErot_CDSD_pc(self, df, calc_Evib_harmonic_anharmonic=False):
        """Calculate Evib & Erot in Lines database:

        - Evib is fetched from Energy Levels database
        - Erot is calculated with Erot = E - Evib

        Note: p, c, j, n is a partition and we just a groupby(these) so all
        poluy, wangu etc. are the same

        Parameters
        ----------

        df: DataFrame
            list of transitions

        Other Parameters
        ----------------

        calc_Evib_harmonic_anharmonic: boolean
            if ``True``, calculate harmonic and anharmonic components of
            vibrational energies (for Treanor distributions)
        """

        # Check inputs
        if calc_Evib_harmonic_anharmonic:
            raise NotImplementedError

        molecule = self.input.molecule
        state = self.input.state  # electronic state
        # TODO: for multi-molecule mode: add loops on molecules and states too
        assert molecule == "CO2"

        # Check energy levels are here
        for iso in self._get_isotope_list(molecule):
            if not iso in self.get_partition_function_molecule(molecule):
                raise AttributeError(
                    "No Partition function calculator defined for isotope {0}".format(
                        iso
                    )
                    + ". You need energies to calculate a non-equilibrium spectrum!"
                    + " Fill the levels parameter in your database definition, "
                    + " with energies of known format: {0}".format(KNOWN_LVLFORMAT)
                    + ". See SpectrumFactory.load_databank() help for more details"
                )

        self.profiler.start("fetch_energy", 2)

        def get_Evib_CDSD_pc_1iso(df, iso):
            """Calculate Evib for a given isotope (energies are specific to a
            given isotope)"""

            # list of energy levels for given isotope
            energies = self.get_energy_levels(molecule, iso, state)

            # only keep vibrational energies
            # see text for how we define vibrational energy
            index = ["p", "c"]
            energies = energies.drop_duplicates(index, inplace=False)
            # (work on a copy)

            # reindexing to get a direct access to level database (instead of using df.v1==v1 syntax)
            energies.set_index(index, inplace=True)

            # Calculate vibrational / rotational levels for all transitions
            def fillEvibu(r):
                """
                Note: we just did a groupby('polyu', 'wangu') so all
                r.poluy, r.wangu etc. are the same in a group r

                Implementation
                ------
                r.polyu.iloc[0],r.wangu.iloc[0],r.ranku.iloc[0]  : [0] because they're
                            all the same
                r.ju.iloc[0]  not necessary (same Tvib) but explicitely mentionning it
                         yields a x20 on performances (60s -> 3s)
                """
                r["Evibu"] = energies.at[(r.polyu.iloc[0], r.wangu.iloc[0]), "Evib"]
                return r

            def fillEvibl(r):
                # Not: p, c, j, n is a partition and we just did a groupby(these)
                # so all r.poluy, r.wangu etc. are the same
                r["Evibl"] = energies.at[(r.polyl.iloc[0], r.wangl.iloc[0]), "Evib"]
                return r

            #        df['Evibu'] = df.groupby(by=['polyu','wangu','ranku']).apply(fillEvibu)
            #        df['Evibl'] = df.groupby('polyl','wangl','rankl').apply(Evibl, axis=1)
            #        %timeit: 43.4s per loop

            # total:  ~ 15s on 460k lines   (probably faster since neq==0.9.20)
            #            try:
            # ~ 6.6 s (probably faster since neq==0.9.20 (radis<1.0)
            df = df.groupby(by=["polyu", "wangu"]).apply(fillEvibu)
            # ~ 6.6 s (probably faster since neq==0.9.20 (radis<1.0)
            df = df.groupby(by=["polyl", "wangl"]).apply(fillEvibl)
            # TODO : use map(dict) version (see _add_EvibErot_CDSD_pcN)

            #            except KeyError:
            #                import traceback
            #                traceback.print_exc()
            #                raise KeyError("{0} -> An error (see above) occured that usually ".format(sys.exc_info()[1]) +
            #                               "happens when the energy level is not referenced in the database. " +
            #                               "Check your partition function calculator, and energies " +
            #                               "for isotope {0} (Factory.parsum_calc['CO2'][{0}]['X'].df)".format(iso))

            # Another version that failed because twice slower than apply() in that case
            # ~ keep it for information
            #        Evibdict = energies.set_index(['p','c','N'])['Evib']
            #        Evibdict = Evibdict.drop_duplicates()
            #        try:
            #            dgb = df.groupby(by=['polyu', 'wangu', 'ranku'])
            #            for (poly, wang, rank), idx in dgb.indices.items(): # ~ 3600 items for 460k lines -> total 15s
            #                Evib = Evibdict[(poly, wang, rank)]             # ~ 7.15 µs
            #                df.loc[idx, 'Evibu'] = Evib                     # ~ 4.38ms
            #
            #            dgb = df.groupby(by=['polyl', 'wangl', 'rankl'])
            #            for (poly, wang, rank), idx in dgb.indices.items(): # ~ 3600 items for 460k lines -> total 15s
            #                Evib = Evibdict[(poly, wang, rank)]             # ~ 7.15 µs
            #                df.loc[idx, 'Evibl'] = Evib                     # ~ 4.38ms

            return df.loc[idx, ["Evibl", "Evibu"]]

        #        df = df.groupby('iso').apply(lambda x: add_Evib_CDSD_pc_1iso(x, x.name))

        df["Evibl"] = np.nan
        df["Evibu"] = np.nan
        for iso, idx in df.groupby("iso").indices.items():
            df.loc[idx, ["Evibl", "Evibu"]] = get_Evib_CDSD_pc_1iso(df.loc[idx], iso)

            if radis.config["DEBUG_MODE"]:
                assert (df.loc[idx, "iso"] == iso).all()

        # Get rotational energy: better recalculate than look up the database
        # (much faster!: perf ~25s -> 765µs)
        df["Erotu"] = df.Eu - df.Evibu
        df["Erotl"] = df.El - df.Evibl

        self.profiler.stop(
            "fetch_energy", f"Fetched energies for all {len(df)} transitions"
        )

        if __debug__:
            self.assert_no_nan(df, "Evibu")
            self.assert_no_nan(df, "Evibl")

        return  # None: Dataframe updated

    def _add_EvibErot_CDSD_pcN(self, df, calc_Evib_harmonic_anharmonic=False):
        """Calculate Evib & Erot in Lines database:

        - Evib is fetched from Energy Levels database
        - Erot is calculated with Erot = E - Evib

        Note: p, c, j, n is a partition and we just a groupby(these) so all
        poluy, wangu etc. are the same

        Parameters
        ----------

        df: DataFrame
            list of transitions

        Other Parameters
        ----------------

        calc_Evib_harmonic_anharmonic: boolean
            if ``True``, calculate harmonic and anharmonic components of
            vibrational energies (for Treanor distributions)
        """
        if __debug__:
            printdbg(
                "called _add_EvibErot_CDSD(calc_Evib_harmonic_anharmonic={0})".format(
                    calc_Evib_harmonic_anharmonic
                )
            )

        if calc_Evib_harmonic_anharmonic:
            raise NotImplementedError

        molecule = self.input.molecule
        state = self.input.state  # electronic state
        # TODO: for multi-molecule mode: add loops on molecules and states too
        assert molecule == "CO2"

        self.profiler.start("fetch_energy_2", 2)

        # Check energy levels are here
        for iso in self._get_isotope_list():
            if not iso in self.get_partition_function_molecule(molecule):
                raise AttributeError(
                    "No Partition function calculator defined for isotope {0}".format(
                        iso
                    )
                    + ". You need energies to calculate a non-equilibrium spectrum!"
                    + " Fill the levels parameter in your database definition, "
                    + " with energies of known format: {0}".format(KNOWN_LVLFORMAT)
                    + ". See SpectrumFactory.load_databank() help for more details"
                )

        def get_Evib_CDSD_pcN_1iso(df, iso):
            """Calculate Evib for a given isotope (energies are specific to a
            given isotope)"""

            # list of energy levels for given isotope
            energies = self.get_energy_levels(molecule, iso, state)

            # only keep vibrational energies
            # see text for how we define vibrational energy
            index = ["p", "c", "N"]
            energies = energies.drop_duplicates(index, inplace=False)
            # (work on a copy)

            # reindexing to get a direct access to level database (instead of using df.v1==v1 syntax)
            energies.set_index(index, inplace=True)
            Evib_dict = dict(list(zip(energies.index, energies.Evib)))

            # Add lower state energy
            df_pcN = df.set_index(["polyl", "wangl", "rankl"])
            df["Evibl"] = df_pcN.index.map(Evib_dict.get).values
            # Add upper state energy
            df_pcN = df.set_index(["polyu", "wangu", "ranku"])
            df["Evibu"] = df_pcN.index.map(Evib_dict.get).values

            return df.loc[:, ["Evibl", "Evibu"]]

        df["Evibl"] = np.nan
        df["Evibu"] = np.nan

        # multiple-isotopes in database
        if "iso" in df:
            for iso, idx in df.groupby("iso").indices.items():
                df.loc[idx, ["Evibl", "Evibu"]] = get_Evib_CDSD_pcN_1iso(
                    df.loc[idx], iso
                )

                if radis.config["DEBUG_MODE"]:
                    assert (df.loc[idx, "iso"] == iso).all()

        else:
            iso = df.attrs["iso"]
            df.loc[:, ["Evibl", "Evibu"]] = get_Evib_CDSD_pcN_1iso(df, iso)

        # Get rotational energy: better recalculate than look up the database
        # (much faster!: perf ~25s -> 765µs)
        df["Erotu"] = df.Eu - df.Evibu
        df["Erotl"] = df.El - df.Evibl

        self.profiler.stop(
            "fetch_energy_2", f"Fetched energies for all {len(df)} transitions"
        )

        if __debug__:
            self.assert_no_nan(df, "Evibu")
            self.assert_no_nan(df, "Evibl")

        return  # None: Dataframe updated

    def _add_EvibErot_CDSD_pcJN(self, df, calc_Evib_harmonic_anharmonic=False):
        """Calculate Evib & Erot in Lines database:

        - Evib is fetched from Energy Levels database
        - Erot is calculated with Erot = E - Evib

        Parameters
        ----------
        df: DataFrame
            list of transitions

        Other Parameters
        ----------------
        calc_Evib_harmonic_anharmonic: boolean
            if ``True``, calculate harmonic and anharmonic components of
            vibrational energies (for Treanor distributions)
        """
        if __debug__:
            printdbg(
                "called _add_EvibErot_CDSD_pcJN(calc_Evib_harmonic_anharmonic={0})".format(
                    calc_Evib_harmonic_anharmonic
                )
            )

        if calc_Evib_harmonic_anharmonic:
            raise NotImplementedError

        molecule = self.input.molecule
        state = self.input.state  # electronic state
        # TODO: for multi-molecule mode: add loops on molecules and states too
        assert molecule == "CO2"

        self.profiler.start("fetch_energy_3", 2)

        # Check energy levels are here
        for iso in self._get_isotope_list(molecule):
            if not iso in self.get_partition_function_molecule(molecule):
                raise AttributeError(
                    "No Partition function calculator defined for isotope {0}".format(
                        iso
                    )
                    + ". You need energies to calculate a non-equilibrium spectrum!"
                    + " Fill the levels parameter in your database definition, "
                    + " with energies of known format: {0}".format(KNOWN_LVLFORMAT)
                    + ". See SpectrumFactory.load_databank() help for more details"
                )

        def get_Evib_CDSD_pcJN_1iso(df, iso):
            """Calculate Evib for a given isotope (energies are specific to a
            given isotope)

            Notes
            -----
            for devs:

            Unlike get_EvibErot_CDSD_pcN_1iso and get_EvibErot_CDSD_pc_1iso,
            no need to use groupby() here, as (per construction) there is only
            one level for a combination of p, c, J, N
            """

            # list of energy levels for given isotope
            energies = self.get_energy_levels(molecule, iso, state)

            # reindexing to get a direct access to level database (instead of using df.v1==v1 syntax)
            index = ["p", "c", "j", "N"]
            #            energies.set_index(index, inplace=True) # cant get it work for some reason
            # ... (it seems groupby().apply() is building some cache variables and
            # ... I cant reset the index of energies properly)
            energies = energies.set_index(index, inplace=False)
            Evib_dict = dict(list(zip(energies.index, energies.Evib)))

            # Add lower state energy
            df_pcJN = df.set_index(["polyl", "wangl", "jl", "rankl"])
            #            for i in df.index:
            #                df.loc[i, 'Evibl'] = energies.at[i, 'Evib']
            # the map below is crazy fast compared to above loop
            df.loc[:, "Evibl"] = df_pcJN.index.map(Evib_dict.get).values
            # Add upper state energy
            df_pcJN = df.set_index(["polyu", "wangu", "ju", "ranku"])
            #            for i in df.index:
            #                df.loc[i, 'Evibu'] = energies.at[i, 'Evib']
            # the map below is crazy fast compared to above loop
            df.loc[:, "Evibu"] = df_pcJN.index.map(Evib_dict.get).values

            return df.loc[:, ["Evibl", "Evibu"]]

        #        df = df.groupby('iso').apply(lambda x: get_Evib_CDSD_pcJN_1iso(x, x.name))

        # slower than the following:
        df["Evibl"] = np.nan
        df["Evibu"] = np.nan

        # multiple-isotopes in database
        if "iso" in df:
            for iso, idx in df.groupby("iso").indices.items():
                df.loc[idx, ["Evibl", "Evibu"]] = get_Evib_CDSD_pcJN_1iso(
                    df.loc[idx], iso
                )

                if radis.config["DEBUG_MODE"]:
                    assert (df.loc[idx, "iso"] == iso).all()

        else:
            iso = df.attrs["iso"]
            df.loc[:, ["Evibl", "Evibu"]] = get_Evib_CDSD_pcJN_1iso(df, iso)

        # Get rotational energy: better recalculate than look up the database
        # (much faster!: perf ~25s -> 765µs)
        df["Erotu"] = df.Eu - df.Evibu
        df["Erotl"] = df.El - df.Evibl

        self.profiler.stop(
            "fetch_energy_3", f"Fetched energies for all {len(df)} transitions"
        )

        if __debug__:
            self.assert_no_nan(df, "Evibu")
            self.assert_no_nan(df, "Evibl")

        return  # None: Dataframe updated

    def _add_Evib123Erot_CDSD_pc(self, df, calc_Evib_harmonic_anharmonic=False):
        """Lookup Evib1, Evib2, Evib3 & Erot for all lines in dataframe.

        Note: p, c, j, n is a partition and we just a groupby(these) so all
        poluy, wangu etc. are the same

        Parameters
        ----------
        df: DataFrame
            list of transitions

        Other Parameters
        ----------------
        calc_Evib_harmonic_anharmonic: boolean
            if ``True``, calculate harmonic and anharmonic components of
            vibrational energies (for Treanor distributions)
        """
        if __debug__:
            printdbg(
                "called _add_Evib123Erot_CDSD_pc(calc_Evib_harmonic_anharmonic={0})".format(
                    calc_Evib_harmonic_anharmonic
                )
            )

        if calc_Evib_harmonic_anharmonic:
            raise NotImplementedError

        molecule = self.input.molecule
        state = self.input.state  # electronic state
        # TODO: for multi-molecule mode: add loops on molecules and states too
        assert molecule == "CO2"

        self.profiler.start(
            "fetch_energy_4",
            2,
            "... Fetching vib123 / rot energies for all {0} transitions".format(
                len(df)
            ),
        )

        # Get Energy database
        if self.parsum_calc == {}:
            raise AttributeError(
                "No Partition function calculator defined in this database"
                + ". You need energies to calculate a non-equilibrium spectrum!"
                + " Fill the levels parameter in your database definition, "
                + " with energies of known format: {0}".format(KNOWN_LVLFORMAT)
                + ". See SpectrumFactory.load_databank() help for more details"
            )

        def get_Evib123_CDSD_pc_1iso(df, iso):
            """Calculate Evib for a given isotope (energies are specific to a
            given isotope)"""
            # TODO: implement with map() instead (much faster!! see get_Evib_CDSD_* )

            energies = self.get_energy_levels(molecule, iso, state)
            # reindexing to get a direct access to level database (instead of using df.v1==v1 syntax)

            # only keep vibrational energies
            # see text for how we define vibrational energy
            index = ["p", "c"]
            energies = energies.drop_duplicates(index, inplace=False)
            # (work on a copy)

            # reindexing to get a direct access to level database (instead of using df.v1==v1 syntax)
            energies.set_index(index, inplace=True)

            # Calculate vibrational / rotational levels for all transitions
            def fillEvib123u(r):
                """
                Note: p, c, j, n is a partition and we just did a groupby(these)
                # so all r.poluy, r.wangu etc. are the same

                Implementation
                ------
                r.polyu.iloc[0],r.wangu.iloc[0],r.ranku.iloc[0]  : [0] because they're
                            all the same
                r.ju.iloc[0]  not necessary (same Tvib) but explicitely mentionning it
                         yields a x20 on performances (60s -> 3s)
                (probably faster since neq==0.9.20) (radis<1.0)
                """
                r["Evib1u"] = energies.at[(r.polyu.iloc[0], r.wangu.iloc[0]), "Evib1"]
                r["Evib2u"] = energies.at[(r.polyu.iloc[0], r.wangu.iloc[0]), "Evib2"]
                r["Evib3u"] = energies.at[(r.polyu.iloc[0], r.wangu.iloc[0]), "Evib3"]
                return r

            def fillEvib123l(r):
                # Not: p, c, j, n is a partition and we just did a groupby(these)
                # so all r.poluy, r.wangu etc. are the same
                r["Evib1l"] = energies.at[(r.polyl.iloc[0], r.wangl.iloc[0]), "Evib1"]
                r["Evib2l"] = energies.at[(r.polyl.iloc[0], r.wangl.iloc[0]), "Evib2"]
                r["Evib3l"] = energies.at[(r.polyl.iloc[0], r.wangl.iloc[0]), "Evib3"]
                return r

            #        df['Evibu'] = df.groupby(by=['polyu','wangu','ranku']).apply(fillEvibu)
            #        df['Evibl'] = df.groupby('polyl','wangl','rankl').apply(Evibl, axis=1)
            #        %timeit: 43.4s per loop

            #            try:  # total:  ~ 15s on 460k lines
            # ~ 6.6 s   (probably faster since neq==0.9.20) (radis<1.0)
            df = df.groupby(by=["polyu", "wangu"]).apply(fillEvib123u)
            # ~ 6.6 s   (probably faster since neq==0.9.20) (radis<1.0)
            df = df.groupby(by=["polyl", "wangl"]).apply(fillEvib123l)
            #            except KeyError:
            #                printr("{0} -> An error (see above) occured that usually ".format(sys.exc_info()[1]) +
            #                       "happens when the energy level is not referenced in the database. " +
            #                       "Check your partition function calculator, and energies " +
            #                       "for isotope {0} (Factory.parsum_calc['CO2'][{0}]['X'].df)".format(iso))
            #                raise

            return df.loc[
                :, ["Evib1l", "Evib2l", "Evib3l", "Evib1u", "Evib2u", "Evib3u"]
            ]

        #        df = df.groupby('iso').apply(lambda x: get_Evib123_CDSD_pc_1iso(x, x.name))

        # Slower than the version below:
        df["Evib1l"] = np.nan
        df["Evib2l"] = np.nan
        df["Evib3l"] = np.nan
        df["Evib1u"] = np.nan
        df["Evib2u"] = np.nan
        df["Evib3u"] = np.nan
        for iso, idx in df.groupby("iso").indices.items():
            df.loc[
                idx, ["Evib1l", "Evib2l", "Evib3l", "Evib1u", "Evib2u", "Evib3u"]
            ] = get_Evib123_CDSD_pc_1iso(df.loc[idx], iso)

        # Add total vibrational energy too (doesnt cost much, and can plot populations in spectrum)
        df["Evibu"] = df.Evib1u + df.Evib2u + df.Evib3u
        df["Evibl"] = df.Evib1l + df.Evib2l + df.Evib3l

        # Get rotational energy: better recalculate than look up the database
        # (much faster! perf:  ~25s -> 765µs)
        df["Erotu"] = df.Eu - df.Evibu
        df["Erotl"] = df.El - df.Evibl

        if __debug__:
            self.assert_no_nan(df, "Evib1u")
            self.assert_no_nan(df, "Evib2u")
            self.assert_no_nan(df, "Evib3u")
            self.assert_no_nan(df, "Evib1l")
            self.assert_no_nan(df, "Evib2l")
            self.assert_no_nan(df, "Evib3l")

        self.profiler.stop("fetch_energy_4", "Fetched energies")
        return  # None: Dataframe updated

    def _add_EvibErot_RADIS_cls1(self, df, calc_Evib_harmonic_anharmonic=False):
        """Fetch Evib & Erot in dataframe for HITRAN class 1 (diatomic)
        molecules.

        Parameters
        ----------
        df: DataFrame
            list of transitions

        Other Parameters
        ----------------
        calc_Evib_harmonic_anharmonic: boolean
            if ``True``, calculate harmonic and anharmonic components of
            vibrational energies (for Treanor distributions)
        """
        if __debug__:
            printdbg(
                "called _add_EvibErot_RADIS_cls1(calc_Evib_harmonic_anharmonic={0})".format(
                    calc_Evib_harmonic_anharmonic
                )
            )

        if calc_Evib_harmonic_anharmonic:
            raise NotImplementedError

        molecule = self.input.molecule
        state = self.input.state  # electronic state
        # TODO: for multi-molecule mode: add loops on molecules and states too

        self.profiler.start("fetch_energy_5", 2)

        # Check energy levels are here
        for iso in self._get_isotope_list():
            if not iso in self.get_partition_function_molecule(molecule):
                raise AttributeError(
                    "No Partition function calculator defined for isotope {0}".format(
                        iso
                    )
                    + ". You need energies to calculate a non-equilibrium spectrum!"
                    + " Fill the levels parameter in your database definition, "
                    + " with energies of known format: {0}".format(KNOWN_LVLFORMAT)
                    + ". See SpectrumFactory.load_databank() help for more details"
                )

        def get_Evib_RADIS_cls1_1iso(df, iso):
            """Calculate Evib & Erot for a given isotope.

            (energies are specific to a given isotope)
            """
            energies = self.get_energy_levels(molecule, iso, state)
            # TODO: for multi-molecule mode: add loops on molecules and states too

            # only keep vibrational energies
            index = ["v"]
            energies = energies.drop_duplicates(index, inplace=False)
            # (work on a copy)

            # reindexing to get a direct access to level database (instead of using df.v1==v1 syntax)
            energies.set_index(index, inplace=True)
            Evib_dict = dict(list(zip(energies.index, energies.Evib)))

            # Add lower state energy
            df_v = df.set_index(["vl"])
            df["Evibl"] = df_v.index.map(Evib_dict.get).values

            # Add upper state energy
            df_v = df.set_index(["vu"])
            df["Evibu"] = df_v.index.map(Evib_dict.get).values

            return df.loc[:, ["Evibl", "Evibu"]]

        df["Evibl"] = np.nan
        df["Evibu"] = np.nan

        # multiple-isotopes in database
        if "iso" in df:
            for iso, idx in df.groupby("iso").indices.items():
                df.loc[idx, ["Evibl", "Evibu"]] = get_Evib_RADIS_cls1_1iso(
                    df.loc[idx], iso
                )

        else:
            iso = df.attrs["iso"]
            df.loc[:, ["Evibl", "Evibu"]] = get_Evib_RADIS_cls1_1iso(df, iso)

        # Get rotational energy: better recalculate than look up the database (much faster!)
        df["Erotu"] = df.Eu - df.Evibu
        df["Erotl"] = df.El - df.Evibl

        assert np.isnan(df.Evibu).sum() == 0
        assert np.isnan(df.Evibl).sum() == 0

        self.profiler.stop(
            "fetch_energy_5", "Fetched energies for all {0} transitions".format(len(df))
        )

        return  # None: Dataframe updated

    def _add_Evib123Erot_RADIS_cls5(self, df):
        """Fetch Evib & Erot in dataframe for HITRAN class 5 (linear triatomic
        with Fermi degeneracy... = CO2! ) molecules

        Parameters
        ----------
        df: DataFrame
            list of transitions

        Other Parameters
        ----------------
        calc_Evib_harmonic_anharmonic: boolean
            if ``True``, calculate harmonic and anharmonic components of
            vibrational energies (for Treanor distributions)

        """
        if __debug__:
            printdbg("called _add_Evib123Erot_RADIS_cls5()")

        molecule = self.input.molecule
        state = self.input.state  # electronic state
        # TODO: for multi-molecule mode: add loops on molecules and states too

        self.profiler.start("fetch_energy_6", 2)

        # Check energy levels are here
        for iso in self._get_isotope_list(molecule):
            if not iso in self.get_partition_function_molecule(molecule):
                raise AttributeError(
                    "No Partition function calculator defined for isotope {0}".format(
                        iso
                    )
                    + ". You need energies to calculate a non-equilibrium spectrum!"
                    + " Fill the levels parameter in your database definition, "
                    + " with energies of known format: {0}".format(KNOWN_LVLFORMAT)
                    + ". See SpectrumFactory.load_databank() help for more details"
                )

        def get_Evib123_RADIS_cls5_1iso(df, iso):
            """Fetch Evib & Erot for a given isotope.

            energies are specific  a given isotope)

            Notes
            -----
            We comb the Line Database with groups of same vibrational level,
            and fetch the corresponding vibrational energy from the Energy Level
            Database.
            """
            # old loop version with groupby.apply() replaced after commit f014007 (15/08/19)

            # Get the Energy Level Database
            energies = self.get_energy_levels(molecule, iso, state)

            # only keep vibrational energies
            # (work on a copy)
            energies = energies.drop_duplicates("viblvl", inplace=False)

            # reindexing to get a direct access to level database (instead of using df.v1==v1 syntax)
            index = ["v1", "v2", "l2", "v3"]
            energies.set_index(index, inplace=True)
            Evib1_dict = dict(list(zip(energies.index, energies.Evib1)))
            Evib2_dict = dict(list(zip(energies.index, energies.Evib2)))
            Evib3_dict = dict(list(zip(energies.index, energies.Evib3)))

            # Add lower state energy
            df_v1v2l2v3 = df.set_index(["v1l", "v2l", "l2l", "v3l"])
            #            for i in df.index:
            #                df.loc[i, 'Evib1l'] = energies.at[i, 'Evib1']
            # the map below is crazy fast compared to above loop
            df["Evib1l"] = df_v1v2l2v3.index.map(Evib1_dict.get).values
            df["Evib2l"] = df_v1v2l2v3.index.map(Evib2_dict.get).values
            df["Evib3l"] = df_v1v2l2v3.index.map(Evib3_dict.get).values
            # TODO @dev # performance: try getting all 3 values at the same time?
            # with:  Evib123_dict = dict(list(zip(energies.index, (energies.Evib1, energies.Evib2, energies.Evib3))))

            # Add upper state energy
            df_v1v2l2v3 = df.set_index(["v1u", "v2u", "l2u", "v3u"])
            df["Evib1u"] = df_v1v2l2v3.index.map(Evib1_dict.get).values
            df["Evib2u"] = df_v1v2l2v3.index.map(Evib2_dict.get).values
            df["Evib3u"] = df_v1v2l2v3.index.map(Evib3_dict.get).values

            return df.loc[
                :, ["Evib1l", "Evib2l", "Evib3l", "Evib1u", "Evib2u", "Evib3u"]
            ]

        #        df = df.groupby('iso').apply(lambda x: get_Evib123_RADIS_cls5_1iso(x, x.name))

        # Slower than the version below:
        df["Evib1l"] = np.nan
        df["Evib2l"] = np.nan
        df["Evib3l"] = np.nan
        df["Evib1u"] = np.nan
        df["Evib2u"] = np.nan
        df["Evib3u"] = np.nan

        # multiple-isotopes in database
        if "iso" in df:
            for iso, idx in df.groupby("iso").indices.items():
                df.loc[
                    idx, ["Evib1l", "Evib2l", "Evib3l", "Evib1u", "Evib2u", "Evib3u"]
                ] = get_Evib123_RADIS_cls5_1iso(df.loc[idx], iso)

        else:
            iso = df.attrs["iso"]
            df.loc[
                :,
                ["Evib1l", "Evib2l", "Evib3l", "Evib1u", "Evib2u", "Evib3u"],
            ] = get_Evib123_RADIS_cls5_1iso(df, iso)

        # Add total vibrational energy too (doesnt cost much, and can plot populations in spectrum)
        df["Evibu"] = df.Evib1u + df.Evib2u + df.Evib3u
        df["Evibl"] = df.Evib1l + df.Evib2l + df.Evib3l

        # Get rotational energy: better recalculate than look up the database (much faster!)
        df["Erotu"] = df.Eu - df.Evibu
        df["Erotl"] = df.El - df.Evibl

        if __debug__:
            self.assert_no_nan(df, "Evib1u")
            self.assert_no_nan(df, "Evib2u")
            self.assert_no_nan(df, "Evib3u")
            self.assert_no_nan(df, "Evib1l")
            self.assert_no_nan(df, "Evib2l")
            self.assert_no_nan(df, "Evib3l")

        self.profiler.stop(
            "fetch_energy_6", "Fetched energies for all {0} transitions".format(len(df))
        )

        return  # None: Dataframe updated

    def _add_Evib123Erot_RADIS_cls5_harmonicanharmonic(self, df):
        """Fetch Evib & Erot in dataframe for HITRAN class 5 (linear triatomic
        with Fermi degeneracy... i.e CO2 ) molecules.

        Parameters
        ----------
        df: DataFrame
            list of transitions
        """
        if __debug__:
            printdbg("called _add_Evib123Erot_RADIS_cls5_harmonicanharmonic()")

        molecule = self.input.molecule
        state = self.input.state  # electronic state
        # TODO: for multi-molecule mode: add loops on molecules and states too

        self.profiler.start("fetch_energy_7", 2)

        # Check energy levels are here
        for iso in self._get_isotope_list(molecule):
            if not iso in self.get_partition_function_molecule(molecule):
                raise AttributeError(
                    "No Partition function calculator defined for isotope {0}".format(
                        iso
                    )
                    + ". You need energies to calculate a non-equilibrium spectrum!"
                    + " Fill the levels parameter in your database definition, "
                    + " with energies of known format: {0}".format(KNOWN_LVLFORMAT)
                    + ". See SpectrumFactory.load_databank() help for more details"
                )

        def get_Evib123_RADIS_cls5_1iso_ah(df, iso):
            """Fetch Evib & Erot for a given isotope (energies are specific to
            a given isotope). Returns harmonic, anharmonic components.

            Notes
            -----

            We comb the Line Database with groups of same vibrational level,
            and fetch the corresponding vibrational energy from the Energy Level
            Database.
            """
            # TODO: implement with map() instead (much faster!! see get_Evib_CDSD_* )

            # Get the Energy Level Database
            energies = self.get_energy_levels(molecule, iso, state)

            # only keep vibrational energies
            energies = energies.drop_duplicates("viblvl", inplace=False)
            # (work on a copy)

            # reindexing to get a direct access to level database (instead of using df.v1==v1 syntax)
            index = ["v1", "v2", "l2", "v3"]
            energies.set_index(index, inplace=True)
            Evib1_h_dict = dict(list(zip(energies.index, energies.Evib1_h)))
            Evib1_a_dict = dict(list(zip(energies.index, energies.Evib1_a)))
            Evib2_h_dict = dict(list(zip(energies.index, energies.Evib2_h)))
            Evib2_a_dict = dict(list(zip(energies.index, energies.Evib2_a)))
            Evib3_h_dict = dict(list(zip(energies.index, energies.Evib3_h)))
            Evib3_a_dict = dict(list(zip(energies.index, energies.Evib3_a)))

            # Add lower state energy
            df_v1v2l2v3 = df.set_index(["v1l", "v2l", "l2l", "v3l"])
            # the map below is crazy fast compared to above loop
            df["Evib1l_h"] = df_v1v2l2v3.index.map(Evib1_h_dict.get).values
            df["Evib1l_a"] = df_v1v2l2v3.index.map(Evib1_a_dict.get).values
            df["Evib2l_h"] = df_v1v2l2v3.index.map(Evib2_h_dict.get).values
            df["Evib2l_a"] = df_v1v2l2v3.index.map(Evib2_a_dict.get).values
            df["Evib3l_h"] = df_v1v2l2v3.index.map(Evib3_h_dict.get).values
            df["Evib3l_a"] = df_v1v2l2v3.index.map(Evib3_a_dict.get).values
            # TODO @dev # performance: try getting all 3 values at the same time?
            # with:  Evib123_dict = dict(list(zip(energies.index, (energies.Evib1, energies.Evib2, energies.Evib3))))

            # Add upper state energy
            df_v1v2l2v3 = df.set_index(["v1u", "v2u", "l2u", "v3u"])
            df["Evib1u_h"] = df_v1v2l2v3.index.map(Evib1_h_dict.get).values
            df["Evib1u_a"] = df_v1v2l2v3.index.map(Evib1_a_dict.get).values
            df["Evib2u_h"] = df_v1v2l2v3.index.map(Evib2_h_dict.get).values
            df["Evib2u_a"] = df_v1v2l2v3.index.map(Evib2_a_dict.get).values
            df["Evib3u_h"] = df_v1v2l2v3.index.map(Evib3_h_dict.get).values
            df["Evib3u_a"] = df_v1v2l2v3.index.map(Evib3_a_dict.get).values

            return df.loc[
                :,
                [
                    "Evib1l_h",
                    "Evib1l_a",
                    "Evib2l_h",
                    "Evib2l_a",
                    "Evib3l_h",
                    "Evib3l_a",
                    "Evib1u_h",
                    "Evib1u_a",
                    "Evib2u_h",
                    "Evib2u_a",
                    "Evib3u_h",
                    "Evib3u_a",
                ],
            ]

        # Slower than the version below:
        df["Evib1l_h"] = np.nan
        df["Evib1l_a"] = np.nan
        df["Evib2l_h"] = np.nan
        df["Evib2l_a"] = np.nan
        df["Evib3l_h"] = np.nan
        df["Evib3l_a"] = np.nan
        df["Evib1u_h"] = np.nan
        df["Evib1u_a"] = np.nan
        df["Evib2u_h"] = np.nan
        df["Evib2u_a"] = np.nan
        df["Evib3u_h"] = np.nan
        df["Evib3u_a"] = np.nan
        for iso, idx in df.groupby("iso").indices.items():
            df.loc[
                idx,
                [
                    "Evib1l_h",
                    "Evib1l_a",
                    "Evib2l_h",
                    "Evib2l_a",
                    "Evib3l_h",
                    "Evib3l_a",
                    "Evib1u_h",
                    "Evib1u_a",
                    "Evib2u_h",
                    "Evib2u_a",
                    "Evib3u_h",
                    "Evib3u_a",
                ],
            ] = get_Evib123_RADIS_cls5_1iso_ah(df.loc[idx], iso)

        # Add total vibrational energy too (doesnt cost much, and can plot populations in spectrum)
        df["Evibu_a"] = df.Evib1u_a + df.Evib2u_a + df.Evib3u_a
        df["Evibu_h"] = df.Evib1u_h + df.Evib2u_h + df.Evib3u_h
        df["Evibl_a"] = df.Evib1l_a + df.Evib2l_a + df.Evib3l_a
        df["Evibl_h"] = df.Evib1l_h + df.Evib2l_h + df.Evib3l_h
        df["Evib1u"] = df.Evib1u_h + df.Evib1u_a
        df["Evib1l"] = df.Evib1l_h + df.Evib1l_a
        df["Evib2u"] = df.Evib2u_h + df.Evib2u_a
        df["Evib2l"] = df.Evib2l_h + df.Evib2l_a
        df["Evib3u"] = df.Evib3u_h + df.Evib3u_a
        df["Evib3l"] = df.Evib3l_h + df.Evib3l_a
        df["Evibu"] = df.Evib1u + df.Evib2u + df.Evib3u
        df["Evibl"] = df.Evib1l + df.Evib2l + df.Evib3l

        # Get rotational energy: better recalculate than look up the database (much faster!)
        df["Erotu"] = df.Eu - df.Evibu
        df["Erotl"] = df.El - df.Evibl

        if __debug__:
            self.assert_no_nan(df, "Evib1u")
            self.assert_no_nan(df, "Evib2u")
            self.assert_no_nan(df, "Evib3u")
            self.assert_no_nan(df, "Evib1l")
            self.assert_no_nan(df, "Evib2l")
            self.assert_no_nan(df, "Evib3l")

        self.profiler.stop(
            "fetch_energy_7", "Fetched energies for all {0} transitions".format(len(df))
        )

        return df

    def _add_ju(self, df):
        """Calculate J'    (upper state)

        Returns
        -------
        None:
            df is updated automatically with column ``'ju'``

        Notes
        -----
        Reminder::

            P branch: J' - J'' = -1
            Q branch: J' - J'' = 0
            R branch: J' - J'' = 1

        P, Q, R are replaced by -1, 0, 1 in the DataFrame to ensure that all
        terms are numeric (improves performances)
        """
        if "branch" not in df:
            raise KeyError(
                f"`branch` not defined in database columns: ({list(df.columns)}). "
                + "You can add it with `load_columns=['branch', ...]` or simply `load_columns='noneq'` / `load_columns='all'` in fetch_databank / load_databank()"
            )

        if df.dtypes["branch"] != np.int64:
            raise DeprecationWarning(
                "For performance purpose, numeric (-1, 0, 1) "
                + "format for (P, Q, R) branch is now required. "
                + "If using cache files, regenerate them?"
            )
        # TODO @dev : switch to CategorialGroup ? (less memory)

        #        df['ju'] = df.jl
        #        df.loc[df.branch==-1,'ju'] -= 1    # branch P
        #        df.loc[df.branch==1,'ju'] += 1     # branch R

        #        # slightly less readable but ~ 20% faster than above:
        dgb = df.groupby("branch")
        df["ju"] = df.jl
        for branch, idx in dgb.indices.items():
            if branch == -1:  # 'P':
                df.loc[idx, "ju"] -= 1
            if branch == 1:  #'R':
                df.loc[idx, "ju"] += 1

        return None

    def _add_Eu(self, df):
        """Calculate upper state energy.

        Returns
        -------
        None:
            df is updated automatically with new column ``'Eu'``
        """

        # Get upper state energy
        df["Eu"] = df.El + df.wav

        return None

    def _calc_noneq_parameters(self, vib_distribution, singleTvibmode):
        """Make sure database has non equilibrium quantities (Evib, Erot, etc.)

        Notes
        -----

        This may be a bottleneck for a first calculation (has to calculate
        the nonequilibrium energies)
        """

        # Checks and loads Energy level database
        if self.misc.load_energies == False:
            self._init_rovibrational_energies(self.levels, self.params.levelsfmt)
            self.misc.load_energies = True

        self.profiler.start("check_non_eq_param", 2)

        df = self.df0
        if len(df) == 0:
            return  # no lines

        # Check spectroscopic parameters required for non-equilibrium (to identify lines)
        for k in ["branch"]:
            if k not in df:
                error_details = ""
                if "globu" in df:
                    error_details = ". However, it looks like `globu` is defined. Maybe HITRAN-like database wasn't fully parsed? See radis.io.hitran.hit2df"
                raise KeyError(
                    f"`{k}` not defined in database ({list(df.columns)}). "
                    + error_details
                    + f"Make sure you properly load parameters required for non-LTE calculations by adding `load_columns=['{k}', ...]` or simply `load_columns='noneq'` in fetch_databank / load_databank()"
                )

        # Make sure database has pre-computed non equilibrium quantities

        # ... Make sure upper J' is calculated  (needed to compute populations)
        if not "ju" in df:
            self._add_ju(df)

        # ... Make sure upper energy level is calculated (needed to compute populations)
        if not "Eu" in df:
            self._add_Eu(df)

        # (Evib, Erot, etc.)
        # This may be a bottleneck for a first calculation (has to calculate
        # the nonequilibrium energies)
        calc_Evib_harmonic_anharmonic = vib_distribution in ["treanor"]
        if calc_Evib_harmonic_anharmonic:
            if singleTvibmode:
                required_columns = ["Evibl_a", "Evibl_h"]
            else:
                required_columns = [
                    "Evib1l_a",
                    "Evib1l_h",
                    "Evib2l_a",
                    "Evib2l_h",
                    "Evib3l_a",
                    "Evib3l_h",
                ]
        else:
            if singleTvibmode:
                required_columns = ["Evibl"]
            else:
                required_columns = ["Evib1l", "Evib2l", "Evib3l"]

        if not all_in(required_columns, df):
            if singleTvibmode:
                self._add_EvibErot(
                    df,
                    calc_Evib_harmonic_anharmonic=calc_Evib_harmonic_anharmonic,
                )
            else:
                self._add_Evib123Erot(
                    df,
                    calc_Evib_harmonic_anharmonic=calc_Evib_harmonic_anharmonic,
                )
            assert all_in(required_columns, df)

        # ... Check no negative energies
        tol = -1e-4  # tolerance for negative energies (in cm-1)
        if not ((df.Erotu > tol).all() and (df.Erotl > tol).all()):
            self.warn(
                "There are negative rotational energies in the database",
                "NegativeEnergiesWarning",
            )

        # ... Make sure degeneracies are calculated
        if not all_in(["gju", "gjl", "gvibu", "gvibl", "gu", "gl"], df):
            self._calc_degeneracies(df)

        if not "Aul" in df:
            self.calc_weighted_trans_moment()
            self.calc_einstein_coefficients()

        self.profiler.stop("check_non_eq_param", "Checked nonequilibrium parameters")

    def _calc_degeneracies(self, df):
        """Calculate vibrational and rotational degeneracies.

        See Also
        --------
        :func:`~radis.db.degeneracies.gs`, :func:`~radis.db.degeneracies.gi`
        """
        from radis.db.degeneracies import gi, gs

        dbformat = self.params.dbformat

        # Rotational + state specific degeneracy (J + isotope dependant)
        df["gju"] = 2 * df.ju + 1
        df["gjl"] = 2 * df.jl + 1

        if "id" in df:
            id_set = df.id.unique()
            if len(id_set) > 1:
                raise NotImplementedError("> 1 molecules in same DataFrame")
            else:
                self.warn(
                    "There shouldn't be a Column 'id' with a unique value",
                    "PerformanceWarning",
                )
            df.attrs["id"] = int(id_set)
            if self.save_memory:
                del df["id"]

        if "iso" in df:  # multiple isotopes in database
            id = df.attrs["id"]
            dgb = df.groupby(by=["iso"])
            for (iso), idx in dgb.indices.items():
                _gs = gs(id, iso)
                if isinstance(_gs, tuple):
                    # Molecules that have alternating degeneracy.
                    if id not in [2]:  # CO2
                        raise NotImplementedError
                    # normally we should find whether the rovibrational level is symmetric
                    # or asymmetric. Here we just assume it's symmetric, because for
                    # symmetric isotopes such as CO2(626), CO2 asymmetric levels
                    # dont exist (gs=0) and they should not be in the line database.
                    _gs = _gs[0]

                dg = df.loc[idx]
                _gi = gi(id, iso)
                df.loc[idx, "grotu"] = dg.gju * _gs * _gi
                df.loc[idx, "grotl"] = dg.gjl * _gs * _gi

                if radis.config["DEBUG_MODE"]:
                    assert (df.loc[idx, "iso"] == iso).all()

        else:
            id = df.attrs["id"]
            isotope = df.attrs["iso"]
            _gs = gs(id, isotope)
            if isinstance(_gs, tuple):
                # Molecules that have alternating degeneracy.
                if id not in [2]:  # CO2
                    raise NotImplementedError
                # normally we should find whether the rovibrational level is symmetric
                # or asymmetric. Here we just assume it's symmetric, because for
                # symmetric isotopes such as CO2(626), CO2 asymmetric levels
                # dont exist (gs=0) and they should not be in the line database.
                _gs = _gs[0]

            dg = df
            _gi = gi(id, isotope)
            df.loc[:, "grotu"] = dg.gju * _gs * _gi
            df.loc[:, "grotl"] = dg.gjl * _gs * _gi

        # %%

        if dbformat in [
            "hitran",
            "hitemp",
            "hitemp-radisdb",
            "cdsd-hitemp",
            "cdsd-4000",
        ]:
            # In HITRAN, AFAIK all molecules have a complete assignment of rovibrational
            # levels hence gvib=1 for all vibrational levels.
            #
            # Complete rovibrational assignment would not be True, for instance, for
            # bending levels of CO2 if all levels are considered degenerated
            # (with a v2+1 degeneracy)
            df["gvibu"] = 1
            df["gvibl"] = 1
        else:
            raise NotImplementedError(
                "vibrational degeneracy assignation for dbformat={0}".format(dbformat)
            )

        # Total
        df["gu"] = df.gvibu * df.grotu
        df["gl"] = df.gvibl * df.grotl

        return None  # dataframe updated directly

    def get_lines_abundance(self, df):
        """Returns the isotopic abundance of each line in `df`

        Parameters
        ----------
        df: dataframe

        Returns
        -------
        float or dict: The abundance of all the isotopes in the dataframe
        """

        molpar = self.molparam

        if "id" in df.columns:
            id_set = df.id.unique()
            if len(id_set) > 1:
                raise NotImplementedError("> 1 molecules in same DataFrame")
            else:
                self.warn(
                    "There shouldn't be a Column 'id' with a unique value",
                    "PerformanceWarning",
                )
                df.attrs["id"] = int(id_set)

        if "iso" in df.columns:
            iso_set = df.iso.unique()
            if len(iso_set) == 1:
                self.warn(
                    "There shouldn't be a Column 'iso' with a unique value",
                    "PerformanceWarning",
                )
                iso = int(iso_set)
                return molpar.get(df.attrs["id"], iso, "abundance")
            else:
                abundance_dict = {}
                for iso in iso_set:
                    abundance_dict[iso] = molpar.get(df.attrs["id"], iso, "abundance")
                return df["iso"].map(abundance_dict)
        else:
            iso = df.attrs["iso"]
            return molpar.get(df.attrs["id"], iso, "abundance")

    def get_molar_mass(self, df):
        """Returns molar mass for all lines of DataFrame ``df``.

        Parameters
        ----------
        df: pd.DataFrame

        Returns
        -------
        The molar mass of all the isotopes in the dataframe


        Notes
        -----

        If molar mass is unknown, it can be temporarily added in the
        :py:attr:`radis.lbl.loader.DatabankLoader._EXTRA_MOLAR_MASS` dictionary
        """
        molpar = self.molparam

        if "id" in df.columns:
            raise NotImplementedError(">1 molecule")
        elif "id" in df.attrs:
            id = df.attrs["id"]
        else:
            # HARDCODED molar mass; for WIP ExoMol implementation, until MolParams
            # is an attribute and can be updated with definitions from ExoMol.
            # https://github.com/radis/radis/issues/321

            # see :py:attr:`radis.lbl.loader.DatabankLoader._EXTRA_MOLAR_MASS`
            try:
                return self._EXTRA_MOLAR_MASS[df.attrs["molecule"]][
                    str(df.attrs["iso"])
                ]
            except KeyError:
                raise NotImplementedError(
                    "Molar mass of {0} (isotope {1}) is unknown.".format(
                        df.attrs["molecule"], df.attrs["iso"]
                    )
                    + " You can manually add it in your radis.json file {'molparams':{'molar_mass':{'molecule':{'ISOTOPE':...}}}}; or in the SpectrumFactory._EXTRA_MOLAR_MASS[molecule][isotope] = M (g/mol) dictionary. Please also report on GitHub so we can update !"
                )

        if "iso" in df.columns:
            iso_set = df.iso.unique()
            molar_mass_dict = {}
            for iso in iso_set:
                molar_mass_dict[iso] = molpar.get(id, iso, "mol_mass")
            molar_mass = df["iso"].map(molar_mass_dict)
        else:
            iso = df.attrs["iso"]
            molar_mass = molpar.get(id, iso, "mol_mass")

        return molar_mass

        #

    def calc_weighted_trans_moment(self):
        """Calculate weighted transition-moment squared :math:`R_s^2` (in ``Debye^2``)

        Returns
        -------
        None:
            ``self.df0`` is updated directly with new column ``Rs2``  .
            R is in ``Debye^2``   (1e-36 ergs.cm3)

        References
        ----------
        Weighted transition-moment squared :math:`R_s^2` from linestrength :math:`S_0`
        at temperature :math:`T_ref`, derived from Eq.(A5) in [Rothman-1998]_

        .. math:
            R_s^2=10^{+36}\\frac{3h c}{8{\\pi}^3} \\frac{1}{n_u} \\frac{1}{\\frac{I_a g_l}{Q_{ref}} \\operatorname{exp}\\left(\\frac{-E_l}{T_{ref}}\\right)} \\frac{1}{1-\\operatorname{exp}\\left(\\frac{-n_u}{T_{ref}}\\right)} S_0
        """

        df = self.df0
        Tref = self.input.Tref

        self.profiler.start("calc_weight_trans", 2)

        # get abundance
        abundance = self.get_lines_abundance(df)
        if not self.molparam.terrestrial_abundances:
            raise NotImplementedError(
                "Formula not corrected for non-terrestrial isotopic abundances"
            )

        gl = df.gl
        El = df.El
        nu = df.wav
        Ia = abundance
        h = h_CGS  # erg.s
        c = c_CGS
        S = (
            df.int
        )  # reference linestrength   ( computed with terrestrial isotopic abundances)

        weighted_trans_moment_sq = (
            (3 * h * c / 8 / pi ** 3)
            / nu
            / (Ia * gl / self.Qgas(df, Tref) * exp(-hc_k * El / Tref))
            / (1 - exp(-hc_k * nu / Tref))
            * 1e36
        ) * S

        df["Rs2"] = weighted_trans_moment_sq

        self.profiler.stop("calc_weight_trans", "calculated weighted transition moment")

        return

    def calc_reference_linestrength(self):
        """Calculate reference linestrength from Einstein coefficients"""

        df = self.df0
        Tref = self.input.Tref

        if self.profiler:
            self.profiler.start("calc_ref_linestrength", 2)

        if "id" in df:
            id_set = df.id.unique()
            if len(id_set) > 1:
                raise NotImplementedError("> 1 molecules in same DataFrame")
            else:
                self.warn(
                    "There shouldn't be a Column 'id' with a unique value",
                    "PerformanceWarning",
                )
            df.attrs["id"] = int(id_set)
        molecule = get_molecule(df.attrs["id"])

        if not "iso" in df:

            # Shortcut if only 1 isotope. We attribute molar_mass & abundance
            # as attributes of the line database, instead of columns. Much
            # faster!

            state = self.input.state
            parsum = self.get_partition_function_calculator(
                molecule, df.attrs["iso"], state
            )  # partition function
            df.attrs["Qref"] = parsum.at(
                Tref, update_populations=False
            )  # stored as attribute, not column
            assert "Qref" not in df.columns
            Qref = df.attrs["Qref"]

        else:
            iso_set = df.iso.unique()
            if len(iso_set) == 1:
                self.warn(
                    "There shouldn't be a Column 'iso' with a unique value",
                    "PerformanceWarning",
                )

            # normal method
            # still much faster than the groupby().apply() method (see radis<=0.9.19)
            # (tested + see https://stackoverflow.com/questions/44954514/efficient-way-to-conditionally-populate-elements-in-a-pandas-groupby-object-pos)

            Qref_dict = {}

            dgb = df.groupby(by=["iso"])
            for (iso), idx in dgb.indices.items():
                state = self.input.state
                parsum = self.get_partition_function_calculator(
                    molecule, iso, state
                )  # partition function
                Qref_dict[iso] = parsum.at(Tref, update_populations=False)
                # ... note: do not update the populations here, so populations in the
                # ... energy level list correspond to the one calculated for T and not Tref

                if radis.config["DEBUG_MODE"]:
                    if "id" in df:
                        assert (df.loc[idx, "id"] == id).all()
                    assert (df.loc[idx, "iso"] == iso).all()

            Qref = df["iso"].map(Qref_dict)
        # NOTE: This S0 is not the same as the one calculated by calc_S0()!!!
        # The difference is that this one is multiplied by the fractional population of transition's levels
        # TO-DO: Resolve this ambiguity
        S0 = linestrength_from_Einstein(
            A=df.A, gu=df.gu, El=df.El, Ia=df.Ia, nu=df.wav, Q=Qref, T=Tref
        )

        if self.profiler:
            self.profiler.stop(
                "calc_ref_linestrength", "Calculated reference linestrength"
            )

        return S0

    def calc_einstein_coefficients(self):
        """Calculate :math:`A_{ul}`, :math:`B_{lu}`, :math:`B_{ul}` Einstein coefficients from weighted
        transition moments squared :math:`R_s^2`.

        Returns
        -------
        None: ``self.df0`` is updated directly with new columns ``Aul``, ``Blu``, ``Bul``

        Notes
        -----
        Einstein A coefficient already in database under df0.A
        Difference between df0.A and df0.Aul < 0.5%

        References
        ----------
        Einstein induced absorption coefficient (in :math:`cm^3/J/s^2`)

        .. math::
            B_{lu}=10^{-36}\\cdot\\frac{8{\\pi}^3}{3h^2} R_s^2 \\cdot 10^{-7}

        Einstein induced emission coefficient (in :math:`cm^3/J/s^2`)

        .. math::
            B_{ul}=10^{-36}\\cdot\\frac{8{\\pi}^3}{3h^2} \\frac{gl}{gu} R_s^2 \\cdot 10^{-7}

        Einstein spontaneous emission coefficient (in :math:`s^{-1}`)

        .. math::
            A_{ul}=10^{-36}\\cdot\\frac{\\frac{64{\\pi}^4}{3h} {\\nu}^3 gl}{gu} R_s^2

        See (Eqs.(A7), (A8), (A9) in [Rothman-1998]_)

        """

        df = self.df0

        try:
            df["Rs2"]
        except KeyError:
            raise KeyError("Weighted transition moment squared not calculated")

        Rs2 = df.Rs2
        gl = df.gl
        gu = df.gu
        nu = df.wav
        h = h_CGS  # erg.s

        # Calculate coefficients
        df["Blu"] = 8 * pi ** 3 / (3 * h ** 2) * Rs2 * 1e-36 * 1e7  # cm3/(J.s^2)
        df["Bul"] = (
            8 * pi ** 3 / (3 * h ** 2) * (gl / gu) * Rs2 * 1e-36 * 1e7
        )  # cm3/(J.s^2)
        df["Aul"] = 64 * pi ** 4 / (3 * h) * nu ** 3 * gl / gu * Rs2 * 1e-36  # s-1

        return None  # dataframe updated directly

    # %% ======================================================================
    # PRIVATE METHODS - APPLY ENVIRONMENT PARAMETERS
    # (all functions that depends upon T or P)
    # (calculates populations, linestrength & radiance, lineshift)
    # (computation: work on df1, called by or after eq_spectrum() )
    # ---------------------------------
    # calc_lineshift
    # calc_linestrength_eq
    # calc_populations_eq
    # calc_populations_noneq
    # calc_linestrength_noneq
    # calc_emission_integral
    # _cutoff_linestrength

    # XXX =====================================================================

    def calc_lineshift(self):
        """Calculate lineshift due to pressure.

        Returns
        -------
        None: ``self.df1`` is updated directly with new column ``shiftwav``

        References
        ----------
        Shifted line center based on pressure shift coefficient :math:`lambda_{shift}`
        and pressure :math:`P`.

        .. math::
            \\omega_{shift}=\\omega_0+\\lambda_{shift} P

        See Eq.(A13) in [Rothman-1998]_
        """

        self.profiler.start("calc_lineshift", 2)

        df = self.df1

        # Calculate
        air_pressure = self.input.pressure_mbar / 1013.25  # convert from mbar to atm

        if "Pshft" in df.columns:
            df["shiftwav"] = df.wav + (df.Pshft * air_pressure)
        else:
            self.warn(
                "Pressure-shift coefficient not given in database: assumed 0 pressure shift",
                "MissingPressureShiftWarning",
            )
            df["shiftwav"] = df.wav

        # Sorted lines is needed for sparse wavenumber range algorithm.
        df.sort_values("shiftwav", inplace=True)

        self.profiler.stop("calc_lineshift", "Calculated lineshift")

        return

    def calc_S0(self):
        """Calculate the unscaled intensity from the tabulated Einstein coefficient.

        Parameters
        ----------
        None

        Returns
        -------
        None: ``self.df0`` is updated directly with new column ``S0``

        References
        ----------

        .. math::
            S_0 = \\frac{I_a g' A_{21}}{8 \\pi c \\nu^2}

        Notes
        -----
        Currently this value is only used in GPU calculations.
        It is one of the columns that is transferred to the GPU
        memory. The idea behind S0 is that it is scaled with all
        variabled that do not change during iterations as to
        minimize calculations.

        Units: cm-1/(molecules/cm-2

        NOTE: S0 is not directly related to S(T) used elsewhere!!!
              (It may even differ in units!!!)

        """

        from radis.phys.constants import c  # m.s-1

        c_cm = c * 100  # cm.s-1

        df0 = self.df0

        if len(df0) == 0:
            return  # no lines

        self.profiler.start("scaled_S0", 2, "... Scaling equilibrium linestrength")

        gp = df0["gp"]
        A = df0["A"]
        wav = df0["wav"]
        Ia = self.get_lines_abundance(df0)

        S0 = Ia * gp * A / (8 * pi * c_cm * wav ** 2)

        df0["S0"] = S0  # [cm-1/(molecules/cm-2)]

        assert "S0" in self.df0

        self.profiler.stop("scaled_S0", "Scaled equilibrium linestrength")

        return

    def _calc_Q(self, molecule, iso, state, T):
        """Get partition function at temperature ``T`` from tabulated values, try with
        calculated partition function (full summation) if Out of Bounds.

        Returns
        -------
        Q: float
            partition functions at temperature ``T``
        """

        try:
            parsum = self.get_partition_function_interpolator(molecule, iso, state)
            Q = parsum.at(T)
        except OutOfBoundError as err:
            # Try to calculate
            try:
                parsum = self.get_partition_function_calculator(molecule, iso, state)
            except KeyError:  # partition function not defined, raise the initial error
                raise err
            else:
                self.warn(
                    "Error with tabulated partition function"
                    + "({0}). Using calculated one instead".format(err.args[0]),
                    "OutOfBoundWarning",
                )
            Q = parsum.at(T)
        return Q

    # Partition functions
    def Qgas(self, df1, Tgas):
        """Calculate partition function Qgas at temperature ``Tgas``, for all lines
        of ``df1``. Returns a single value if all lines have the same Qgas value,
        or a column if they are different

        Parameters
        ----------
        Tgas: float (K)
            gas temperature

        Returns
        -------
        float or dict: Returns Qgas as a dictionary with isotope values as its keys

        See Also
        --------
        :py:meth:`~radis.lbl.base.BaseFactory.Qgas_Qref_ratio`
        """

        if "id" in df1.columns:
            id_set = df1.id.unique()
            if len(id_set) > 1:
                raise NotImplementedError(">1 molecule.")
            else:
                self.warn(
                    "There shouldn't be a Column 'id' with a unique value",
                    "PerformanceWarning",
                )
                df1.attrs["id"] = int(id_set)

        if "molecule" in df1.attrs:
            molecule = df1.attrs["molecule"]  # used for ExoMol, which has no HITRAN-id
        else:
            molecule = get_molecule(df1.attrs["id"])
        state = self.input.state

        if "iso" in df1:
            iso_set = df1.iso.unique()
            if len(iso_set) == 1:
                self.warn(
                    "There shouldn't be a Column 'iso' with a unique value",
                    "PerformanceWarning",
                )
                iso = int(iso_set)
                Q = self._calc_Q(molecule, iso, state, Tgas)
                df1.attrs["Q"] = Q
                return Q
            else:
                Qgas_dict = {}
                for iso in iso_set:
                    Qgas_dict[iso] = self._calc_Q(molecule, iso, state, Tgas)
                return df1["iso"].map(Qgas_dict)

        else:  # "iso" not in df:
            iso = df1.attrs["iso"]
            Q = self._calc_Q(molecule, iso, state, Tgas)
            df1.attrs["Q"] = Q
            return Q

    def Qref_Qgas_ratio(self, df1, Tgas, Tref):
        """Calculate Qref/Qgas at temperature ``Tgas``, ``Tref``, for all lines
        of ``df1``. Returns a single value if all lines have the same Qref/Qgas ratio,
        or a column if they are different

        See Also
        --------
        :py:meth:`~radis.lbl.base.BaseFactory.Qgas`
        """

        if "id" in df1.columns:
            id_set = df1.id.unique()
            if len(id_set) > 1:
                raise NotImplementedError(">1 molecule.")
            else:
                self.warn(
                    "There shouldn't be a Column 'id' with a unique value",
                    "PerformanceWarning",
                )
                df1.attrs["id"] = int(id_set)

        if "molecule" in df1.attrs:
            molecule = df1.attrs["molecule"]  # used for ExoMol, which has no HITRAN-id
        else:
            molecule = get_molecule(df1.attrs["id"])
        state = self.input.state

        if "iso" in df1:
            iso_set = df1.iso.unique()
            if len(iso_set) == 1:
                self.warn(
                    "There shouldn't be a Column 'iso' with a unique value",
                    "PerformanceWarning",
                )
                iso = int(iso_set)
                Qgas = self._calc_Q(molecule, iso, state, Tgas)
                Qref = self._calc_Q(molecule, iso, state, Tref)
                df1.attrs["Qgas"] = Qgas
                df1.attrs["Qref"] = Qref
                Qref_Qgas = Qref / Qgas

            else:
                Qref_Qgas_ratio = {}
                for iso in iso_set:
                    Qgas = self._calc_Q(molecule, iso, state, Tgas)
                    Qref = self._calc_Q(molecule, iso, state, Tref)
                    Qref_Qgas_ratio[iso] = Qref / Qgas
                Qref_Qgas = df1["iso"].map(Qref_Qgas_ratio)

        else:
            iso = df1.attrs["iso"]
            Qgas = self._calc_Q(molecule, iso, state, Tgas)
            Qref = self._calc_Q(molecule, iso, state, Tref)
            df1.attrs["Qgas"] = Qgas
            df1.attrs["Qref"] = Qref
            Qref_Qgas = Qref / Qgas
        return Qref_Qgas

    def calc_linestrength_eq(self, Tgas):
        """Calculate linestrength at temperature Tgas correcting the database
        linestrength tabulated at temperature :math:`T_{ref}`.

        Parameters
        ----------
        Tgas: float (K)
            gas temperature

        Returns
        -------
        None: ``self.df1`` is updated directly with new column ``S``

        References
        ----------

        .. math::
            S(T) = S_0 \\frac{Q_{ref}}{Q_{gas}} \\operatorname{exp}\\left(-E_l \\left(\\frac{1}{T_{gas}}-\\frac{1}{T_{ref}}\\right)\\right) \\frac{1-\\operatorname{exp}\\left(\\frac{-\\omega_0}{Tgas}\\right)}{1-\\operatorname{exp}\\left(\\frac{-\\omega_0}{T_{ref}}\\right)}

        See Eq.(A11) in [Rothman-1998]_

        Notes
        -----
        Internals:

        (some more informations about what this function does)

        Starts with df1 which is still a copy of df0 loaded by
        :meth:`~radis.lbl.loader.DatabankLoader.load_databank`
        Updates linestrength in df1. Cutoff criteria is applied afterwards.

        .. minigallery:: radis.lbl.base.BaseFactory.calc_linestrength_eq
            :add-heading:

        See Also
        --------
        :py:func:`~radis.lbl.base.linestrength_from_Einstein`
        """

        Tref = self.input.Tref
        df1 = self.df1

        if len(df1) == 0:
            return  # no lines

        self.profiler.start(
            "scaled_eq_linestrength", 2, "... Scaling equilibrium linestrength"
        )

        # %% Calculate line strength at desired temperature
        # -------------------------------------------------

        if self.molparam.terrestrial_abundances:

            # This calculation is based on equation (A11) in Rothman 1998: "JQSRT, vol.
            # 60, No. 5, pp. 665-710"

            # correct for Partition Function
            df1["S"] = (
                df1.int
                * self.Qref_Qgas_ratio(df1, Tgas, Tref)
                *
                # ratio of Boltzman populations
                exp(-hc_k * df1.El * (1 / Tgas - 1 / Tref))
                *
                # effect of stimulated emission
                (1 - exp(-hc_k * df1.wav / Tgas))
                / (1 - exp(-hc_k * df1.wav / Tref))
            )  # [cm-1/(molecules/cm-2)]

        else:
            # An alternative strategy is to calculate the linestrength from the
            # Einstein A coefficient and the populations (see Klarenaar 2017 Eqn. 12)

            if not "gu" in df1:
                if not "ju" in df1:
                    self._add_ju(df1)
                self._calc_degeneracies(df1)

            Ia = self.get_lines_abundance(df1)
            df1["S"] = linestrength_from_Einstein(
                df1.A, df1.gu, df1.El, Ia, df1.wav, self.Qgas(df1, Tgas), Tgas
            )

        assert "S" in self.df1

        self.profiler.stop("scaled_eq_linestrength", "Scaled equilibrium linestrength")

        return

    # %%
    def calc_populations_eq(self, Tgas):
        """Calculate upper state population for all active transitions in
        equilibrium case (only used in total power calculation)

        Parameters
        ----------
        Tgas: float (K)
            temperature

        Returns
        -------
        None:
            `nu` is stored in self.df1

        Notes
        -----
        Isotopes: these populations are not corrected for the isotopic abundance,
        i.e, abundance has to be accounted for if used for emission density
        calculations (based on Einstein A coefficient), but not for linestrengths
        (that include the abundance dependency already)

        References
        ----------
        Population of upper state follows a Boltzmann distribution:

        .. math::
            n_u = g_u \\frac{\\operatorname{exp}\\left(\\frac{-E_u}{T_{gas}}\\right)}{Q_{gas}}

        See Also
        --------
        :meth:`~radis.lbl.base.BaseFactory.calc_populations_noneq`,
        :meth:`~radis.lbl.base.BaseFactory._calc_populations_noneq_multiTvib`,
        :meth:`~radis.levels.partfunc.RovibPartitionFunction.at`
        """

        df1 = self.df1

        self.profiler.start("calc_eq_population", 2)

        # Calculate degeneracies
        # ----------------------------------------------------------------------

        if not "ju" in df1:
            self._add_ju(df1)

        if not "gu" in df1:
            self._calc_degeneracies(df1)

        # ... Make sure upper energy level is calculated (needed to compute populations)
        if not "Eu" in df1:
            self._add_Eu(df1)

        # Calculate population
        # ----------------------------------------------------------------------
        # equilibrium: Boltzmann in any case
        df1["nu"] = (
            df1.gu.values * exp(-hc_k * df1.Eu.values / Tgas) / self.Qgas(df1, Tgas)
        )

        assert "nu" in self.df1

        self.profiler.stop("calc_eq_population", "Calculated equilibrium populations")

        return

    def Qneq(self, df, Tvib, Trot, vib_distribution, rot_distribution, overpopulation):
        """Nonequilibrium partition function

        Returns
        -------
        column or float: depending if there are many isotopes or one"""
        self.profiler.start("part_function", 3)

        if "id" in df:
            id_set = df.id.unique()
            if len(id_set) > 1:
                raise NotImplementedError("> 1 molecules in same DataFrame")
            else:
                self.warn(
                    "There shouldn't be a Column 'id' with a unique value",
                    "PerformanceWarning",
                )
                df.attrs["id"] = int(id_set)
        molecule = get_molecule(df.attrs["id"])
        state = self.input.state

        if "iso" in df:
            Q_dict = {}
            iso_set = df.iso.unique()
            if len(iso_set) == 1:
                self.warn(
                    "There shouldn't be a Column 'iso' with a unique value",
                    "PerformanceWarning",
                )
            for iso in iso_set:
                parsum = self.get_partition_function_calculator(molecule, iso, state)
                if is_float(Tvib):
                    Q_dict[iso] = parsum.at_noneq(
                        Tvib,
                        Trot,
                        vib_distribution=vib_distribution,
                        rot_distribution=rot_distribution,
                        overpopulation=overpopulation,
                        update_populations=self.misc.export_populations,
                    )
                else:
                    Q_dict[iso] = parsum.at_noneq_3Tvib(
                        Tvib,
                        Trot,
                        vib_distribution=vib_distribution,
                        rot_distribution=rot_distribution,
                        overpopulation=overpopulation,
                        update_populations=self.misc.export_populations,
                    )
            Q = df["iso"].map(Q_dict)

        else:  # "iso" not in df:
            iso = df.attrs["iso"]
            parsum = self.get_partition_function_calculator(molecule, iso, state)
            if is_float(Tvib):
                Q = parsum.at_noneq(
                    Tvib,
                    Trot,
                    vib_distribution=vib_distribution,
                    rot_distribution=rot_distribution,
                    overpopulation=overpopulation,
                    update_populations=self.misc.export_populations,
                )
            else:
                Q = parsum.at_noneq_3Tvib(
                    Tvib,
                    Trot,
                    vib_distribution=vib_distribution,
                    rot_distribution=rot_distribution,
                    overpopulation=overpopulation,
                    update_populations=self.misc.export_populations,
                )
            df.attrs["Q"] = Q
        self.profiler.stop("part_function", "partition functions")
        return Q

    def Qneq_Qvib_Qrotu_Qrotl(
        self, df, Tvib, Trot, vib_distribution, rot_distribution, overpopulation
    ):
        """Nonequilibrium partition function; with the detail of
        vibrational partition function and rotational partition functions"""
        # TODO @ dev : implement the map(dict) approach to fill Q Qvib Qrotu Qrotl
        # note : Qrot already use a map(dict) so we need a map with 2 keys. Is it worth it?
        self.profiler.start("part_function", 3)

        if "id" in df:
            id_set = df.id.unique()
            if len(id_set) > 1:
                raise NotImplementedError("> 1 molecules in same DataFrame")
            else:
                self.warn(
                    "There shouldn't be a Column 'id' with a unique value",
                    "PerformanceWarning",
                )
                df.attrs["id"] = int(id_set)
        molecule = get_molecule(df.attrs["id"])
        state = self.input.state

        if "iso" in df:  #  multiple isotopes
            iso_set = df.iso.unique()
            if len(iso_set) == 1:
                self.warn(
                    "There shouldn't be a Column 'iso' with a unique value",
                    "PerformanceWarning",
                )

            dgb = df.groupby(by=["iso"])
            for (iso), idx in dgb.indices.items():

                # Get partition function for all lines
                parsum = self.get_partition_function_calculator(molecule, iso, state)

                if is_float(Tvib):
                    Q, Qvib, dfQrot = parsum.at_noneq(
                        Tvib,
                        Trot,
                        vib_distribution=vib_distribution,
                        rot_distribution=rot_distribution,
                        overpopulation=overpopulation,
                        returnQvibQrot=True,
                        update_populations=self.misc.export_populations,
                    )
                else:
                    raise NotImplementedError(
                        "Cannot return detail of Qvib, Qrot for 3-Tvib mode"
                    )

                # ... make sure PartitionFunction above is calculated with the same
                # ... temperatures, rovibrational distributions and overpopulations
                # ... as the populations of active levels (somewhere below)
                df.at[idx, "Qvib"] = Qvib
                df.at[idx, "Q"] = Q

                # reindexing to get a direct access to Qrot database
                # create the lookup dictionary
                # dfQrot index is already 'viblvl'
                dfQrot_dict = dict(list(zip(dfQrot.index, dfQrot.Qrot)))

                dg = df.loc[idx]

                # Add lower state Qrot
                dg_sorted = dg.set_index(["viblvl_l"], inplace=False)
                df.loc[idx, "Qrotl"] = dg_sorted.index.map(dfQrot_dict.get).values
                # Add upper state energy
                dg_sorted = dg.set_index(["viblvl_u"], inplace=False)
                df.loc[idx, "Qrotu"] = dg_sorted.index.map(dfQrot_dict.get).values

                if radis.config["DEBUG_MODE"]:
                    assert (df.loc[idx, "iso"] == iso).all()

            Q, Qvib, Qrotu, Qrotl = df.Q, df.Qvib, df.Qrotu, df.Qrotl

        else:
            iso = df.attrs["iso"]

            parsum = self.get_partition_function_calculator(molecule, iso, state)

            if is_float(Tvib):
                Q, Qvib, dfQrot = parsum.at_noneq(
                    Tvib,
                    Trot,
                    vib_distribution=vib_distribution,
                    rot_distribution=rot_distribution,
                    overpopulation=overpopulation,
                    returnQvibQrot=True,
                    update_populations=self.misc.export_populations,
                )
            else:
                raise NotImplementedError(
                    "Cannot return detail of Qvib, Qrot for 3-Tvib mode"
                )

            # ... make sure PartitionFunction above is calculated with the same
            # ... temperatures, rovibrational distributions and overpopulations
            # ... as the populations of active levels (somewhere below)
            df.attrs["Qvib"] = Qvib
            df.attrs["Q"] = Q
            assert "Qvib" not in df.columns
            assert "Q" not in df.columns

            # reindexing to get a direct access to Qrot database
            # create the lookup dictionary
            # dfQrot index is already 'viblvl'
            dfQrot_dict = dict(list(zip(dfQrot.index, dfQrot.Qrot)))

            dg = df.loc[:]

            # Add lower state Qrot
            dg_sorted = dg.set_index(["viblvl_l"], inplace=False)
            df.loc[:, "Qrotl"] = dg_sorted.index.map(dfQrot_dict.get).values
            # Add upper state energy
            dg_sorted = dg.set_index(["viblvl_u"], inplace=False)
            df.loc[:, "Qrotu"] = dg_sorted.index.map(dfQrot_dict.get).values

            Q, Qvib, Qrotu, Qrotl = Q, Qvib, df.Qrotu, df.Qrotl

        self.profiler.stop("part_function", "partition functions")
        return Q, Qvib, Qrotu, Qrotl

    # %%
    def calc_populations_noneq(
        self,
        Tvib,
        Trot,
        vib_distribution="boltzmann",
        rot_distribution="boltzmann",
        overpopulation=None,
    ):
        """Calculate upper and lower state population for all active
        transitions, as well as all levels (through
        :meth:`~radis.levels.partfunc.RovibPartitionFunction.at_noneq`)

        Parameters
        ----------
        Tvib, Trot: float (K)
            temperatures
        vib_distribution: ``'boltzmann'``, ``'treanor'``
            vibrational level distribution
        rot_distribution: ``'boltzmann'``
            rotational level distribution
        overpopulation: dict, or ``None``
            dictionary of overpopulation factors for vibrational levels

        Returns
        -------
        None: `nu`, `nl`, `nu_vib`, `nl_vib` are stored in self.df1

        Notes
        -----
        Isotopic abundance:

        Note that these populations are not corrected for the isotopic abundance,
        i.e, abundance has to be accounted for if used for emission density
        calculations (based on Einstein A coefficient), but not for linestrengths
        (that include the abundance dependency already)

        All populations:

        This method calculates populations of emitting and absorbing levels.
        Populations of all levels (even the one not active on the spectral
        range considered) are calculated during the Partition function calculation.
        See: :meth:`~radis.levels.partfunc.RovibPartitionFunction.at_noneq`

        References
        ----------
        Boltzmann vibrational distributions

        .. math::

            n_{vib}=\\frac{g_{vib}}{Q_{vib}} \\operatorname{exp}\\left(\\frac{-E_{vib}}{T_{vib}}\\right)

        or Treanor vibrational distributions

        .. math::

            n_{vib}=\\frac{g_{vib}}{Qvib} \\operatorname{exp}\\left(-\\left(\\frac{E_{vib,harm}}{T_{vib}}+\\frac{E_{vib,anharm}}{T_{rot}}\\right)\\right)

        Overpopulation of vibrational levels

        .. math::

            n_{vib}=\\alpha n_{vib}

        Boltzmann rotational distributions

        .. math::

            n_{rot}=\\frac{g_{rot}}{Q_{rot}} \\operatorname{exp}\\left(\\frac{-E_{rot}}{T_{rot}}\\right)

        Final rovibrational population of one level

        .. math::

            n=n_{vib} n_{rot} \\frac{Q_{rot} Q_{vib}}{Q}

        See Also
        --------
        :meth:`~radis.lbl.base.BaseFactory.calc_populations_eq`,
        :meth:`~radis.lbl.base.BaseFactory._calc_populations_noneq_multiTvib`,
        :meth:`~radis.levels.partfunc.RovibPartitionFunction.at_noneq`
        """

        # Check inputs
        if overpopulation is None:
            overpopulation = {}

        df = self.df1

        self.profiler.start("calc_noneq_population", 2)

        if len(df) == 0:
            return  # no lines in database, no need to go further

        # Get vibrational levels for both upper and lower states
        if not ("viblvl_u" in df and not "viblvl_l" in df):
            from radis.lbl.bands import add_bands

            add_bands(
                df,
                dbformat=self.params.dbformat,
                lvlformat=self.params.levelsfmt,
                verbose=self.verbose,
            )
            assert "viblvl_u" in df
            assert "viblvl_l" in df

        # %%

        #  Derive populations
        if not self.misc.export_rovib_fraction:
            if overpopulation != {}:
                raise NotImplementedError(
                    "Overpopulation not implemented in multi-Tvib mode"
                )
            # ... vibrational distribution
            if vib_distribution == "boltzmann":
                df["nu_vib_x_Qvib"] = df.gvibu * exp(-hc_k * df.Evibu / Tvib)
                df["nl_vib_x_Qvib"] = df.gvibl * exp(-hc_k * df.Evibl / Tvib)
            elif vib_distribution == "treanor":
                raise NotImplementedError("TO DO!")  #!!!TODO
            else:
                raise ValueError(
                    "Unknown vibrational distribution: {0}".format(vib_distribution)
                )

            # ... Rotational distributions
            if rot_distribution == "boltzmann":
                df["nu_rot_x_Qrot"] = df.grotu * exp(-df.Erotu * hc_k / Trot)
                df["nl_rot_x_Qrot"] = df.grotl * exp(-df.Erotl * hc_k / Trot)
            else:
                raise ValueError(
                    "Unknown rotational distribution: {0}".format(rot_distribution)
                )

            # ... Partition functions
            Qneq = self.Qneq(
                df,
                Tvib,
                Trot,
                vib_distribution=vib_distribution,
                rot_distribution=rot_distribution,
                overpopulation=overpopulation,
            )

            # ... Total
            df["nu"] = df.nu_vib_x_Qvib * df.nu_rot_x_Qrot / Qneq
            df["nl"] = df.nl_vib_x_Qvib * df.nl_rot_x_Qrot / Qneq

        else:  # self.misc.export_rovib_fraction:
            # ... Partition functions
            Q, Qvib, Qrotu, Qrotl = self.Qneq_Qvib_Qrotu_Qrotl(
                df,
                Tvib,
                Trot,
                vib_distribution=vib_distribution,
                rot_distribution=rot_distribution,
                overpopulation=overpopulation,
            )

            # ... vibrational distribution
            if vib_distribution == "boltzmann":
                # equation generated with @pytexit.py2tex > see docstrings.
                df["nu_vib"] = df.gvibu / Qvib * exp(-hc_k * df.Evibu / Tvib)
                df["nl_vib"] = df.gvibl / Qvib * exp(-hc_k * df.Evibl / Tvib)
            elif vib_distribution == "treanor":
                df["nu_vib"] = (
                    df.gvibu
                    / Qvib
                    * exp(-hc_k * (df.Evibu_h / Tvib + df.Evibu_a / Trot))
                )
                df["nl_vib"] = (
                    df.gvibl
                    / Qvib
                    * exp(-hc_k * (df.Evibl_h / Tvib + df.Evibl_a / Trot))
                )
            else:
                raise ValueError(
                    "Unknown vibrational distribution: {0}".format(vib_distribution)
                )

            # ... Add vibrational-specific overpopulation factors
            if overpopulation != {}:
                for viblvl, ov in overpopulation.items():
                    if ov != 1:
                        df.loc[df.viblvl_u == viblvl, "nu_vib"] *= ov
                        df.loc[df.viblvl_l == viblvl, "nl_vib"] *= ov

            # ... Rotational distributions
            if rot_distribution == "boltzmann":
                df["nu_rot"] = df.grotu / df.Qrotu * exp(-df.Erotu * hc_k / Trot)
                df["nl_rot"] = df.grotl / df.Qrotl * exp(-df.Erotl * hc_k / Trot)
            else:
                raise ValueError(
                    "Unknown rotational distribution: {0}".format(rot_distribution)
                )

            # ... Total
            df["nu"] = df.nu_vib * df.nu_rot * (Qrotu * Qvib / Q)
            df["nl"] = df.nl_vib * df.nl_rot * (Qrotl * Qvib / Q)

        if __debug__:
            assert "nu" in self.df1
            assert "nl" in self.df1
            self.assert_no_nan(self.df1, "nu")
            self.assert_no_nan(self.df1, "nl")

        self.profiler.stop(
            "calc_noneq_population", "Calculated nonequilibrium populations"
        )

        return

    # %%
    def _calc_populations_noneq_multiTvib(
        self,
        Tvib,
        Trot,
        vib_distribution="boltzmann",
        rot_distribution="boltzmann",
        overpopulation=None,
    ):
        """Calculate upper and lower state population for all active
        transitions, as well as all levels (through
        :meth:`~radis.levels.partfunc.RovibPartitionFunction.at_noneq`)

        Parameters
        ----------
        Tvib, Trot: float (K)
            temperatures
        vib_distribution: ``'boltzmann'``, ``'treanor'``
            vibrational level distribution
        rot_distribution: ``'boltzmann'``
            rotational level distribution
        overpopulation: dict, or ``None``
            dictionary of overpopulation factors for vibrational levels

        Notes
        -----
        Isotopic abundance:

        Note that these populations are not corrected for the isotopic abundance,
        i.e, abundance has to be accounted for if used for emission density
        calculations (based on Einstein A coefficient), but not for linestrengths
        (that include the abundance dependency already)

        All populations:

        This method calculates populations of emitting and absorbing levels.
        Populations of all levels (even the one not active on the spectral
        range considered) are calculated during the Partition function calculation.
        :meth:`~radis.levels.partfunc.RovibPartitionFunction.at_noneq`

        Todo someday:

        - so far it's a 3 Tvib model, hardcoded. make it a N-vibrational model,
          with lists / dictionary?

        See Also
        --------
        :meth:`~radis.lbl.base.BaseFactory.calc_populations_eq`,
        :meth:`~radis.lbl.base.BaseFactory.calc_populations_noneq`,
        :meth:`~radis.levels.partfunc.RovibPartitionFunction.at_noneq_3Tvib`
        """

        # Check inputs
        if overpopulation is None:
            raise NotImplementedError(
                "Overpopulation not implemented in multi-Tvib mode"
            )
        Tvib1, Tvib2, Tvib3 = Tvib

        df = self.df1

        self.profiler.start("calc_noneq_population_multiTvib", 2)

        if len(df) == 0:
            return  # no lines in database, no need to go further

        # partition function
        # ... unlike the (tabulated) equilibrium case, here we recalculate it from
        # scratch

        #  Derive populations
        # ... vibrational distribution
        if vib_distribution == "boltzmann":
            nu_vib1Qvib1 = df.gvibu * exp(-hc_k * df.Evib1u / Tvib1)
            nl_vib1Qvib1 = df.gvibl * exp(-hc_k * df.Evib1l / Tvib1)
            nu_vib2Qvib2 = df.gvibu * exp(-hc_k * df.Evib2u / Tvib2)
            nl_vib2Qvib2 = df.gvibl * exp(-hc_k * df.Evib2l / Tvib2)
            nu_vib3Qvib3 = df.gvibu * exp(-hc_k * df.Evib3u / Tvib3)
            nl_vib3Qvib3 = df.gvibl * exp(-hc_k * df.Evib3l / Tvib3)
        elif vib_distribution == "treanor":
            # fmt: off
            nu_vib1Qvib1 = df.gvibu * exp(-hc_k * (df.Evib1u_h / Tvib1 + df.Evib1u_a / Trot))
            nl_vib1Qvib1 = df.gvibl * exp(-hc_k * (df.Evib1l_h / Tvib1 + df.Evib1l_a / Trot))
            nu_vib2Qvib2 = df.gvibu * exp(-hc_k * (df.Evib2u_h / Tvib2 + df.Evib2u_a / Trot))
            nl_vib2Qvib2 = df.gvibl * exp(-hc_k * (df.Evib2l_h / Tvib2 + df.Evib2l_a / Trot))
            nu_vib3Qvib3 = df.gvibu * exp(-hc_k * (df.Evib3u_h / Tvib3 + df.Evib3u_a / Trot))
            nl_vib3Qvib3 = df.gvibl * exp(-hc_k * (df.Evib3l_h / Tvib3 + df.Evib3l_a / Trot))
            # fmt: on
        else:
            raise ValueError(
                "Unknown vibrational distribution: {0}".format(vib_distribution)
            )

        if overpopulation != {}:
            raise NotImplementedError(overpopulation)
        # Not Implemented:
        #        if overpopulation != {}:
        #            if not ('viblvl_u' in df and not 'viblvl_l' in df):
        #                from radis.lbl.bands import add_bands
        #                df = add_bands(df, dbformat=self.params.dbformat, verbose=self.verbose)
        #                assert 'viblvl_u' in df
        #                assert 'viblvl_l' in df
        #
        #            for viblvl, ov in overpopulation.items():
        #                if ov != 1:
        #                    df.loc[df.viblvl_u==viblvl, 'nu_vib'] *= ov
        #                    df.loc[df.viblvl_l==viblvl, 'nl_vib'] *= ov

        # ... Rotational distributions
        # that would require Qrot, which we dont have (NotImplemented
        # for 3 temperatures). Let's just get the total

        # ... Total
        if rot_distribution == "boltzmann":
            # ... total
            Qneq = self.Qneq(
                df,
                Tvib,
                Trot,
                vib_distribution=vib_distribution,
                rot_distribution=rot_distribution,
                overpopulation=overpopulation,
            )

            df["nu"] = (
                nu_vib1Qvib1
                * nu_vib2Qvib2
                * nu_vib3Qvib3
                * df.grotu
                * exp(-df.Erotu * hc_k / Trot)
                / Qneq
            )
            df["nl"] = (
                nl_vib1Qvib1
                * nl_vib2Qvib2
                * nl_vib3Qvib3
                * df.grotl
                * exp(-df.Erotl * hc_k / Trot)
                / Qneq
            )

        else:
            raise ValueError(
                "Unknown rotational distribution: {0}".format(rot_distribution)
            )

        assert "nu" in self.df1
        assert "nl" in self.df1

        self.profiler.stop(
            "calc_noneq_population_multiTvib",
            "Calculated nonequilibrium populations (multiTvib)",
        )

        return

    def get_lines(self):
        """Return lines if self.misc.export_lines is True, else get None."""

        if self.misc.export_lines:
            return self.df1
        else:
            return None

    # %% Get populations
    def get_populations(self, levels="vib"):
        """For all molecules / isotopes / electronic states, lookup energy
        levels as calculated in partition function calculators, and (if
        calculated) populations, and returns as a dictionary.

        Parameters
        ----------

        levels: ``'vib'``, ``'rovib'``, list of these, or ``None``
            what levels to get. Note that ``'rovib'`` can yield large Spectrum objects.

        Returns
        -------

        pops: dict
            Structure::

                {molecule: {isotope: {electronic_state: {'vib': pandas Dataframe,    # (copy of) vib levels
                                                         'rovib': pandas Dataframe,  # (copy of) rovib levels
                                                         'Ia': float    # isotopic abundance
                                                         }}}}

        See Also
        --------

        :meth:`~radis.lbl.base.BaseFactory.get_energy_levels`
        """

        # Check input
        if levels is None or levels is False:
            return {}
        if isinstance(levels, str):
            levels = [levels]
        for l in levels:
            EXPECTED = ["vib", "rovib"]
            if l not in EXPECTED:
                raise ValueError(
                    "Unexpect type of levels to return {0}. Expected one of {1}".format(
                        l, EXPECTED
                    )
                )

        # To get isotopic abundance
        # placeholder # TODO: replace with attributes of Isotope>ElectronicState objects
        molpar = self.molparam

        pops = {}
        # Loop over molecules, isotopes, electronic states
        for molecule in [self.input.molecule]:
            pops[molecule] = {}

            for isotope in self._get_isotope_list(molecule):
                pops[molecule][isotope] = {}

                id = get_molecule_identifier(molecule)
                # fetch all table directly
                params = molpar.df.loc[(id, isotope)]
                # placeholder # TODO: replace with attributes of Isotope>ElectronicState objects
                Ia = params.abundance

                for electronic_state in list(self.input.state):
                    pops[molecule][isotope][electronic_state] = {}

                    # Add vib or rovib levels
                    energies = self.get_energy_levels(
                        molecule, isotope, electronic_state
                    )
                    for level_type in levels:

                        if level_type == "vib":
                            assert "viblvl" in list(energies.keys())
                            # only get one entry per vibrational level
                            pop = energies.drop_duplicates("viblvl")  # is a copy
                            # remove unecessary keys (all rotational specific)
                            for k in ["E", "j", "gj", "Erot", "grot", "n"]:
                                try:
                                    del pop[k]
                                except KeyError:
                                    pass

                        elif level_type == "rovib":
                            pop = energies.copy()  # is a copy

                        else:
                            raise ValueError(
                                "Unknown level type: {0}".format(level_type)
                            )

                        # Store
                        pops[molecule][isotope][electronic_state][level_type] = pop

                    # Add extra information (isotope - electronic state specific)
                    pops[molecule][isotope][electronic_state]["Ia"] = Ia

        # Note: all dataframes should be copies, else it gets dangerous if
        # exported in a Spectrum but still connected to the Factory
        return pops

    def calc_linestrength_noneq(self):
        """Calculate linestrengths at non-LTE

        Parameters
        ----------
        Pre-requisite:

            lower state population `nl` has already been calculated by
            :meth:`~radis.lbl.base.BaseFactory.calc_populations_noneq`


        Returns
        -------
        None
            Linestrength `S` added in self.df


        Notes
        -----

        Internals:

        (some more informations about what this function does)

        Starts with df1 which is was a copy of df0 loaded by load_databank(),
        with non-equilibrium quantities added and populations already calculated.
        Updates linestrength in df1. Cutoff criteria is applied afterwards.

        See Also
        --------

        :py:meth:`~radis.lbl.base.BaseFactory.calc_populations_noneq`,
        :py:meth:`~radis.lbl.base.BaseFactory.calc_emission_integral`,
        :py:func:`~radis.lbl.base.linestrength_from_Einstein`

        """

        df = self.df1
        Tref = self.input.Tref

        if len(df) == 0:
            return  # no lines in database, no need to go further

        self.profiler.start("scaled_non_eq_linestrength", 2)

        try:
            df["nl"]
        except KeyError:
            raise KeyError("Calculate populations first")

        # Correct linestrength

        # ... populations without abundance dependance (already in linestrength)
        nu = df.nu  # Note: populations are not corrected for abundance
        nl = df.nl
        # ... remove Qref, nref, etc.
        # start from for tabulated linestrength
        line_strength = df.int.copy()  # TODO: savememory; replace the "int" column
        if not self.molparam.terrestrial_abundances:
            raise NotImplementedError(
                "Formula not corrected for non-terrestrial isotopic abundances"
            )

        # ... correct for lower state population
        line_strength /= df.gl * exp(-hc_k * df.El / Tref) / self.Qgas(df, Tref)
        line_strength *= nl

        # ... correct effect of stimulated emission
        line_strength /= 1 - exp(-hc_k * df.wav / Tref)
        line_strength *= 1 - df.gl / df.gu * nu / nl
        df["S"] = line_strength

        self.profiler.stop(
            "scaled_non_eq_linestrength", "scaled nonequilibrium linestrength"
        )

        return  # df1 automatically updated

    # %%
    def calc_emission_integral(self):
        r"""Calculate Emission Integral.

        .. math::
            Ei=\frac{n_u A_{ul}}{4} \pi \Delta E_{ul}

        Emission Integral is a non usual quantity introduced here to have an
        equivalent of Linestrength in emission calculation

        Returns
        -------
        None
            Emission integral `Ei` added in self.df

        Notes
        -----

        emission_integral: (mW/sr)
            emission integral is defined as::

            Ei = n_u * A_ul / 4π * DeltaE_ul
                  :      :     :      :
                (#/#)  (s-1) (sr)    (mJ)

            So that the radiance ϵ is:

                ϵ(λ)   =      Ei  *   Phi(λ)  * ntot   * path_length
                 :             :         :       :          :
            mW/cm2/sr/nm    (mW/sr)   (1/nm)   (cm-3)      (cm)

        See Also
        --------
        :py:meth:`~radis.lbl.base.BaseFactory.calc_linestrength_noneq`
        """

        df = self.df1

        self.profiler.start("calc_emission_integral", 2)

        if len(df) == 0:
            return  # no lines in database, no need to go further

        try:
            df["nu"]
        except KeyError:
            raise KeyError("Calculate populations first")

        # Calculation

        # adim. (#/#) (multiplied by n_tot later)
        n_u = df["nu"]
        # correct for abundance
        n_ua = n_u * self.get_lines_abundance(df)

        A_ul = df["Aul"]  # (s-1)

        DeltaE = cm2J(df.wav)  # (cm-1) -> (J)
        Ei = n_ua * A_ul / 4 / pi * DeltaE  # (W/sr)

        Ei *= 1e3  # (W/sr) -> (mW/sr)
        df["Ei"] = Ei

        self.profiler.stop("calc_emission_integral", "calculated emission integral")

        return

    # %%
    def _cutoff_linestrength(self, cutoff=None):
        """Discard linestrengths that are lower that this, to reduce
        calculation times. Set the number of lines cut in
        ``self._Nlines_cutoff``

        Parameters
        ----------

        cutoff: float (unit of linestrength:  cm-1/(molecule.cm-2))
            discard linestrengths that are lower that this, to reduce calculation
            times. If 0, no cutoff. Default 0

        Notes
        -----

        # TODO:

        turn linestrength cutoff criteria in 'auto' mode that adjusts linestrength
        calculations based an error percentage criteria
        """

        # Update defaults
        if cutoff is not None:
            self.params.cutoff = cutoff

        # Load variables
        cutoff = self.params.cutoff
        verbose = self.verbose
        df = self.df1

        if len(df) == 0:  # no lines
            self._Nlines_cutoff = None
            return

        if cutoff <= 0:
            self._Nlines_cutoff = 0
            return  # dont update self.df1

        self.profiler.start("applied_linestrength_cutoff", 2)

        # Cutoff:
        b = df.S <= cutoff
        Nlines_cutoff = b.sum()

        # Estimate time gained
        # TODO: Add a better formula to estimate time gained during broadening process
        """ Previous method used was:
            expected_broadening_time_gain = (
                self._broadening_time_ruleofthumb * Nlines_cutoff * len(self.wbroad_centered)
            )
        """

        # Estimate error being made:
        if self.warnings["LinestrengthCutoffWarning"] != "ignore":

            error = df.S[b].sum() / df.S.sum() * 100

            if verbose >= 2:
                print(
                    "Discarded {0:.2f}% of lines (linestrength<{1}cm-1/(#.cm-2))".format(
                        Nlines_cutoff / len(df.S) * 100, cutoff
                    )
                    + " Estimated error: {0:.2f}%".format(error)
                )
            if error > self.misc.warning_linestrength_cutoff:
                self.warn(
                    "Estimated error after discarding lines is large: {0:.2f}%".format(
                        error
                    )
                    + ". Consider reducing cutoff",
                    "LinestrengthCutoffWarning",
                )

        try:
            assert sum(~b) > 0
        except AssertionError as err:
            self.plot_linestrength_hist(cutoff=cutoff)
            raise AssertionError(
                f"All lines discarded! Please increase cutoff (currently : {cutoff:.1e}). "
                + "In your case: (min,max,mean)=({0:.2e},{1:.2e},{2:.2e}".format(
                    df.S.min(), df.S.max(), df.S.mean()
                )
                + "cm-1/(#.cm-2)). See histogram"
            ) from err

        # update df1:
        self.df1 = pd.DataFrame(df[~b])
        #        df.drop(b.index, inplace=True)   # performance: was not faster
        # ... @dev performance: quite long to select here, but I couldn't find a faster
        # ... alternative
        # TODO: remove useless columns in df1 to save memory
        # Note @EP : with Vaex; the selection should be updated here

        # Ensures abundance, molar mass and partition functions are transfered
        # (needed if they are attributes and not isotopes)
        transfer_metadata(df, self.df1, [k for k in df_metadata if k in df.attrs])
        # assert len(self.df1.attrs) > 0

        # Store number of lines cut (for information)
        self._Nlines_cutoff = Nlines_cutoff

        self.profiler.stop("applied_linestrength_cutoff", "Applied linestrength cutoff")

        return

    # %% ======================================================================
    # PRIVATE METHODS - UTILS
    # (cleaning)
    # ---------------------------------
    # _reinitialize_factory
    # _check_inputs
    # plot_populations()

    # XXX =====================================================================

    def _reinitialize(self):
        """Reinitialize Factory before a new spectrum is calculated. It does:

        - create new line Dataframe ``df1`` that will be scaled later with new populations
        - clean some objects if needed to save memory
        - delete populations from RovibrationalPartitionFunction objects

        If in save_memory mode, removes the line database (``self.df0``). This
        function is called after the scaled line database (``self.df1``) has been
        created.
        It saves a lot of memory but prevents the user from calculating a new
        spectrum without reloading the database.

        Returns
        -------

        None:
            but creates ``self.df1`` from ``self.df0``
        """
        if __debug__:
            printdbg("called ._clean_factory()")

        self.profiler.start("reinitialize", 2)

        keep_initial_database = not self.save_memory

        self.profiler.start("copy_database", 3)

        if keep_initial_database:

            # Create new line Dataframe
            # ... Operate on a duplicate dataframe to make it possible to do different
            # ... runs without reloading database
            self.df1 = self.df0.copy()

            # abundance and molar_mass should have been copied even if they are attributes
            # (only 1 molecule, 1 isotope) and not a column (line specific) in the database
            # @dev: this brings a lot of performance improvement, but sometimes fail.
            # | here we ensure that the DataFrame has the values:
            transfer_metadata(
                self.df0, self.df1, [k for k in df_metadata if hasattr(self.df0, k)]
            )

        else:
            self.df1 = self.df0  # self.df0 will be deleted
            del self.df0  # delete attribute name
        self.profiler.stop("copy_database", "Copying database")

        self.profiler.start("memory_usage_warning", 3)
        # Check memory size
        try:
            # Retrieving total user RAM
            mem = virtual_memory()
            mem = mem.total  # total physical memory available
            # Checking if object type column exists
            if "O" in self.df1.dtypes.unique():
                limit = mem / 25  # 4% of user RAM
                self.warn(
                    "'object' type column found in database, calculations and "
                    + "memory usage would be faster with a numeric type. Possible "
                    + "solution is to not use 'save_memory' and convert the columns to dtype.",
                    "PerformanceWarning",
                )
            else:
                limit = mem / 2  # 50 % of user RAM

            # Note: the difference between deep=True and deep=False is around 4 times

            df_size = self.df1.memory_usage(deep=False).sum()

            if df_size > limit:
                self.warn(
                    "Line database is large: {0:.0f} Mb".format(df_size * 1e-6)
                    + ". Consider using save_memory "
                    + "option, if you don't need to reuse this factory to calculate new spectra",
                    "MemoryUsageWarning",
                )
        except ValueError:  # had some unexplained ValueError: __sizeof__() should return >= 0
            pass

        self.profiler.stop("memory_usage_warning", "Check Memory usage of database")

        self.profiler.start("reset_population", 3)
        # Reset populations from RovibrationalPartitionFunctions objects
        molecule = self.input.molecule
        state = self.input.state
        for isotope in self._get_isotope_list(molecule):
            # ... Get partition function calculator
            try:
                parsum = self.get_partition_function_calculator(
                    molecule, isotope, state
                )
            except KeyError:
                # Partition function calculator not defined but maybe it wont
                # be needed (ex: only equilibrium calculations). Continue
                if __debug__:
                    printdbg(
                        "parsum[{0}][{1}][{2}]".format(molecule, isotope, state)
                        + " not defined."
                    )
            else:
                # ... Reset it
                parsum.reset_populations()
        self.profiler.stop("reset_population", "Reset populations")

        self.profiler.stop("reinitialize", "Reinitialize database")

    def _check_inputs(self, mole_fraction, Tmax):
        """Check spectrum inputs, add warnings if suspicious values.

        Also check that line databases look appropriate for the temperature Tmax
        considered

        Parameters
        ----------

        Tmax: float
            Tgas at equilibrium, or max(Tgas, Tvib, Trot) at nonequilibrium
        """

        # Check mole fractions
        if mole_fraction is None:
            raise ValueError("Set mole_fraction")
        if mole_fraction > 1:
            self.warn("mole_fraction is > 1", "InputConditionsWarning")

        # Check temperature range
        if Tmax > 700 and self.params.dbformat in ["hitran"]:
            self.warn(
                "HITRAN is valid for low temperatures (typically < 700 K). "
                + "For higher temperatures you may need HITEMP or CDSD. See the "
                + "'databank=' parameter",
                "HighTemperatureWarning",
            )

        # Also check some computation parameters:

        # Checks there if there is change in wstep value if initial wstep != "auto" (could happen if users modified the _wstep value directly)
        if self._wstep != "auto":
            assert self._wstep == self.params.wstep
        if self._sparse_ldm != "auto":
            assert self._sparse_ldm == self.params.sparse_ldm

        # Checks there if there is change in truncation value
        # (except in the case where truncation is None, where we set it to be the full range)
        if self.params.truncation is not None:
            assert self.truncation == self.params.truncation

        # Check neighbour lines wasn't changed since first initialisation
        # (can create problems if database is not reloaded
        if self._neighbour_lines != self.params.neighbour_lines:
            raise AssertionError(
                f"neighbour_lines value changed from {self._neighbour_lines} to "
                + f"{self.params.neighbour_lines}. Did you reset it manually ? This is currently forbidden as new "
                + "lines won't be retrieved from the database"
            )
            # note @dev:  could be implemented; i.e. send `neighbour_lines` to load_databank instead of SpectrumFactory initialisation

    def plot_populations(self, what="vib", isotope=None, nfig=None):
        """Plot populations currently calculated in factory.

        Plot populations of all levels that participate in the partition function.
        Output is different from the
        Spectrum :py:meth:`~radis.spectrum.spectrum.Spectrum.plot_populations` method,
        where only the levels that directly contribute to the spectrum are shown

        Note: only valid after calculating non_eq spectrum as it uses the
        partition function calculator object

        Parameters
        ----------

        what: 'vib', 'rovib'
            what levels to plot

        isotope: int, or ``None``
            which isotope to plot. If ``None`` and if there are more than one isotope,
            raises an error.

        Other Parameters
        ----------------

        nfig: int, or str
            on which Figure to plot. Default ``None``
        """

        import matplotlib.pyplot as plt

        # Check inputs
        assert what in ["vib", "rovib"]
        if isotope is None:
            if type(self.input.isotope) is int:  # only one isotope. Use it.
                isotope = self.input.isotope
            else:
                raise ValueError("isotope number is needed")

        # Get levels
        molecule = self.input.molecule
        state = self.input.state
        levels = self.get_partition_function_calculator(molecule, isotope, state).df
        if what == "vib":
            E, n, g = levels["Evib"], levels["nvib"], levels["gvib"]
        elif what == "rovib":
            E, n = levels["E"], levels["n"]
            if "g" in list(levels.keys()):
                g = levels["g"]
            elif all_in(["gvib", "grot"], list(levels.keys())):
                g = levels["gvib"] * levels["grot"]
            else:
                raise ValueError(
                    "either g, or gvib+grot must be defined to "
                    + "calculate total degeneracy. Got: {0}".format(list(levels.keys()))
                )

        # Plot
        set_style()
        plt.figure(num=nfig)
        plt.plot(E, n / g, "ok")
        plt.xlabel("Energy (cm-1)")
        plt.ylabel("Population (n / g)")
        plt.yscale("log")
        fix_style()


def get_waverange(
    wmin=None,
    wmax=None,
    wunit=Default("cm-1"),
    wavenum_min=None,
    wavenum_max=None,
    wavelength_min=None,
    wavelength_max=None,
    medium="air",
):
    """Returns wavenumber based on whatever input was given: either ν_min,
    ν_max directly, or λ_min, λ_max  in the given propagation ``medium``.

    Parameters
    ----------
    medium: ``'air'``, ``'vacuum'``
        propagation medium
    wmin, wmax: float, or `~astropy.units.quantity.Quantity` or ``None``
        hybrid parameters that can serve as both wavenumbers or wavelength depending on the unit accompanying them.
        If unitless, wunit is assumed as the accompanying unit.
    wunit: string
        The unit accompanying wmin and wmax. Cannot be passed without passing values for wmin and wmax.
        Default: cm-1
    wavenum_min, wavenum_max: float, or `~astropy.units.quantity.Quantity` or ``None``
        wavenumbers
    wavelength_min, wavelength_max: float, or `~astropy.units.quantity.Quantity` or ``None``
        wavelengths in given ``medium``
    Returns
    -------
    wavenum_min, wavenum_max,: float
        wavenumbers
    """

    # Checking consistency of all input variables

    w_present = wmin is not None and wmax is not None
    wavenum_present = wavenum_min is not None and wavenum_max is not None
    wavelength_present = wavelength_min is not None and wavelength_max is not None
    if w_present + wavenum_present + wavelength_present != 1:
        _msg = "Got : "
        _msg += f" wmin={wmin}, wmax={wmax}" if w_present else ""
        _msg += (
            f" wavenum_min={wavenum_min}, wavenum_max={wavenum_max}"
            if wavenum_present
            else ""
        )
        _msg += (
            f" wavelength_min={wavelength_min}, wavelength_max={wavelength_max}"
            if wavelength_present
            else ""
        )
        raise ValueError(
            "Please pass exactly one range for wavenumber/wavelength input. {}".format(
                _msg
            )
            + "\nChoose either `wmin=..., wmax=...` (with astropy.units), "
            "`wavenum_min=..., wavenum_max=...` (in cm-1), "
            "or `wavelength_min=..., wavelength_max=...` (in nm)"
            "We recommend to use units. Example: \n\n  import astropy.units as u\n  calc_spectrum(wmin=2000 / u.cm, wmax=2300 / u.cm, ..."
        )

    if not isinstance(wunit, Default):
        if not u.Unit(wunit).is_equivalent(u.m) and not u.Unit(wunit).is_equivalent(
            1 / u.m
        ):
            raise ValueError("Wunit dimensions should be either [length] or 1/[length]")

    if wavelength_min is not None or wavelength_max is not None:
        assert wavelength_min is not None and wavelength_max is not None
        assert wavelength_min < wavelength_max
        if not isinstance(wunit, Default):
            raise ValueError("Please use wmin/wmax when passing wunit")
        if isinstance(wavelength_min, u.Quantity):
            assert isinstance(wavelength_min, u.Quantity) and isinstance(
                wavelength_max, u.Quantity
            )
            assert wavelength_min.unit.is_equivalent(u.m)
            assert wavelength_max.unit.is_equivalent(u.m)

    if wavenum_min is not None or wavenum_max is not None:
        assert wavenum_min is not None and wavenum_max is not None
        assert wavenum_min < wavenum_max
        if not isinstance(wunit, Default):
            raise ValueError("Please use wmin/wmax when passing wunit")
        if isinstance(wavenum_min, u.Quantity) or isinstance(wavenum_max, u.Quantity):
            assert isinstance(wavenum_min, u.Quantity) and isinstance(
                wavenum_max, u.Quantity
            )
            assert wavenum_min.unit.is_equivalent(1 / u.cm)
            assert wavenum_max.unit.is_equivalent(1 / u.cm)

    if isinstance(wmin, u.Quantity) or isinstance(wmax, u.Quantity):
        assert wmin is not None and wmax is not None
        assert isinstance(wmin, u.Quantity) and isinstance(wmax, u.Quantity)
        assert wmin.unit.is_equivalent(u.m) or wmin.unit.is_equivalent(1 / u.m)
        assert wmax.unit.is_equivalent(u.m) or wmax.unit.is_equivalent(1 / u.m)
        assert wmin.unit.is_equivalent(wmax.unit)
        if not isinstance(wunit, Default):
            if not wmin.unit.is_equivalent(u.Unit(wunit)):
                raise ValueError("Conflicting units passed for wmin/wmax and wunit")

    # Conversion to base units
    if wmin is not None:
        if isinstance(wmin, u.Quantity) or isinstance(wmax, u.Quantity):
            if wmin.unit.is_equivalent(u.m):
                wavelength_min = wmin
                wavelength_max = wmax
            else:
                wavenum_min = wmin
                wavenum_max = wmax
        else:
            if isinstance(wunit, Default):
                wavenum_min = wmin * u.Unit(wunit.value)
                wavenum_max = wmax * u.Unit(wunit.value)
            else:
                if u.Unit(wunit).is_equivalent(u.m):
                    wavelength_min = wmin * u.Unit(wunit)
                    wavelength_max = wmax * u.Unit(wunit)
                else:
                    wavenum_min = wmin * u.Unit(wunit)
                    wavenum_max = wmax * u.Unit(wunit)

    # We now have wavenum_min/max, or wavelength_min/max defined. Let's convert these to cm-1 (warning: propagating medium is required if we start from wavelengths!)
    if wavenum_min is not None or wavenum_max is not None:
        wavenum_min = convert_and_strip_units(wavenum_min, 1 / u.cm)
        wavenum_max = convert_and_strip_units(wavenum_max, 1 / u.cm)

    elif wavelength_min is not None or wavelength_max is not None:
        if not medium in ["air", "vacuum"]:
            raise NotImplementedError(
                "Unknown propagating medium: {0}. Choose `'air'` or `'vacuum'` to get the correct wavelengths".format(
                    medium
                )
            )
        wavelength_min = convert_and_strip_units(wavelength_min, u.nm)
        wavelength_max = convert_and_strip_units(wavelength_max, u.nm)
        if medium == "air":
            wavenum_min = nm_air2cm(wavelength_max)
            wavenum_max = nm_air2cm(wavelength_min)
        else:  # medium == 'vacuum':
            wavenum_min = nm2cm(wavelength_max)
            wavenum_max = nm2cm(wavelength_min)

    return wavenum_min, wavenum_max


def linestrength_from_Einstein(
    A,
    gu,
    El,
    Ia,
    nu,
    Q,
    T,
):
    r"""Calculate linestrength at temperature ``T`` from Einstein coefficients.

    Parameters
    ----------
    A : float, s-1
        Einstein emission coefficients
    gu : int
        upper state degeneracy
    El : float, cm-1
        lower state energy
    Ia : float
        isotope abundance
    nu : cm-1
        transition wavenumber
    Q : float
        partition function at temperature ``T``
    T : float
        temperature

    Returns
    -------
    S : float
        linestrength at temperature ``T``.

    References
    ----------

    .. math::
        S(T) =\frac{1}{8\pi c_{CGS} {n_u}^2} A \frac{I_a g_u \operatorname{exp}\left(\frac{-c_2 E_l}{T}\right)}{Q(T)} \left(1-\operatorname{exp}\left(\frac{-c_2 n_u}{T}\right)\right)

    Combine Eq.(A.5), (A.9) in [Rothman-1998]_

    See Also
    --------
    :py:meth:`~radis.lbl.base.BaseFactory.calc_linestrength_eq`

    """

    return (
        (1 / (8 * np.pi * c_CGS * nu ** 2))
        * A
        * ((Ia * gu * np.exp(-hc_k * El / T)) / Q)
        * (1 - np.exp(-hc_k * nu / T))
    )


if __name__ == "__main__":
    from radis.test.lbl.test_base import _run_testcases

    _run_testcases()
