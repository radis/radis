# -*- coding: utf-8 -*-
"""Variants of default Partition function tabulators and calculators, based on
CDSD-4000 format, specifically for CO2.

Routine Listing
---------------

Partition functions:

- :class:`~radis.levels.partfunc_cdsd.PartFuncCO2_CDSDcalc`
- :class:`~radis.levels.partfunc_cdsd.PartFuncCO2_CDSDtab`

Which inherit from:

- :class:`~radis.levels.partfunc.RovibParFuncCalculator`
- :class:`~radis.levels.partfunc.RovibParFuncTabulator`


-------------------------------------------------------------------------------
"""


import time
from os.path import exists, getmtime
from warnings import warn

import pandas as pd
from scipy.interpolate import splev, splrep

import radis
from radis.db.molecules import ElectronicState
from radis.io.cache_files import filter_metadata, load_h5_cache_file, save_to_hdf
from radis.lbl.labels import (
    vib_lvl_name_cdsd_p,
    vib_lvl_name_cdsd_pc,
    vib_lvl_name_cdsd_pcJN,
    vib_lvl_name_cdsd_pcN,
)
from radis.levels.partfunc import RovibParFuncCalculator, RovibParFuncTabulator
from radis.misc.warning import OutOfBoundError
from radis.test.utils import getTestFile

# %% Variants of Tabulated partition functions (interpolate)


class PartFuncCO2_CDSDtab(RovibParFuncTabulator):
    """Return partition function of CO2 using a spline interpolation of
    tabulated values used in [CDSD-4000]_

    Parameters
    ----------

    T: temperature (K)
        gas temperature (equilibrium)

    Notes
    -----

    Partition function calculated in CDSD by direct summation (Jmax=300)

    The CDSD ``partition_functions.txt`` can be downloaded from the
    [CDSD-4000]_ FTP : ftp://ftp.iao.ru/pub/CDSD-4000/

    See Also
    --------

    :py:class:`~radis.levels.partfunc_cdsd.PartFuncCO2_CDSDcalc`
    """

    def __init__(self, isotope, database):
        """Get partition function for one isotope only.

        (note that this forces to reload the file once per isotope, but
        at least we have a clean layout with one object per isotope)
        """

        Itable = {1: "12C16O2", 2: "13C16O2", 3: "16O12C18O", 4: "16O12C17O"}

        # Check input
        if not isotope in Itable:
            raise KeyError(
                "Isotope {0} not defined in CDSD tabulated ".format(isotope)
                + "partition functions. Only the following are: {0}".format(
                    list(Itable.keys())
                )
            )

        # Read partition function tabulated data
        parsum = pd.read_csv(database, comment="#", delim_whitespace=True)
        if not "T(K)" in list(parsum.keys()):
            raise KeyError("Missing columns ({0}) in {1}".format("T(K)", database))

        # Define spline interpolation
        isoname = Itable[isotope]
        tck = splrep(parsum["T(K)"], parsum[isoname])

        self.molecule = 2  # 'CO2'
        self.iso = isotope

        self.tck = tck  # dictionary
        self.Tmin = parsum["T(K)"].min()
        self.Tmax = parsum["T(K)"].max()

    def _inrange(self, T):
        """Allow for 5% extrapolation (ex: 296K / 300K) )"""
        return (self.Tmin * 0.95 <= T) and (self.Tmax * 1.05 >= T)

    def _at(self, T):
        """Get partition function at temperature T.

        Called by :meth:`radis.levels.partfunc.RovibParFuncTabulator.at`
        """
        if not self._inrange(T):
            raise OutOfBoundError(
                "Temperature: {0} is out of bounds {1}-{2}".format(
                    T, self.Tmin, self.Tmax
                )
            )

        return splev(T, self.tck)


# %% Calculated partition functions (either from energy levels, or ab initio)


class PartFuncCO2_CDSDcalc(RovibParFuncCalculator):
    """Calculate Partition Function from energy levels (and maybe export a
    tabulated database).

    warning::
        ZPE (zero point energy) must be the same in the Line Database and the
        Energy levels database. See the
        :ref:`Configuration file <label_lbl_config_file>`.

    Parameters
    ----------
    energy_levels: filename
        path to energy levels (to calculate Partition Function) for ``isotope``
    isotope: int
        which isotope we're dealing with. Default ``1``. In the current implementation
        only isotope 1 and 2 are defined.
    levelsfmt: ``'cdsd-p'``, ``'cdsd-pc'``, ``'cdsd-pcN'``, ``'cdsd-hamil'``, or ``None``
        the format of the Energy Database, and in particular how ``Evib`` and ``Erot``
        have been calculated. A vibrational level in the CDSD (p,c,J,N) nomenclature
        can be defined for levels that share a same (p), (p,c) or (p,c,N), where
        ``p`` is the polyad number, ``c`` is the Wang symmetry number, and ``N``
        is the ranking index of a (p,c,J) group. Default ``'cdsd-pc'``.

        If ``None``, dont label the levels. Wont be able to use the EnergyDatabase to fetch
        vibrational energies for lines, however it can still be used to
        calculate Partition functions independently from a Spectrum calculation

    Other Parameters
    ----------------
    use_cached: ``True``, ``False``, or ``'regen'``, ``'force'``
        if ``True``, use (and generate if doesnt exist) a ``.h5`` file.
        If ``'regen'``, regenerate cache file. If ``'force'``, raise an error
        if file doesnt exist. Default ``True``
    mode: 'full summation', 'tabulation'
        calculation mode. ``'tabulation'`` is much faster but not all possible
        distributions are implemented. Default ``'full-summation'``
    use_json: boolean
        deprecated. Better use h5 now.

    Notes
    -----

    Database format:

    Taskhkun database updated with ranking number (n) & total rank (N) of
    block, Evib and Erot (cm-1)  and jref

    For nonequilibrium, different strategies exist so as how to assign rotational and vibrational
    energies in a CDSD database. See the E. Pannier "Limits on CO2 Nonequilibrium
    model" article for a discussion on that.

    Example of table format::

        # calculated rovibrational energy levels of 12C16O2
        # =================================================
        # S.Tashkun, Zuev Institute of Atmospheric Optics, Tomsk, Russia
        # date: 17.03.2017
        #
        # zero point energy ZPE (cm-1) =  2531.828
        #
        # p = 2v1 + v2 + 3v3 - polyad number
        # j - rotational quantum number
        # c - Wang symmetry (1-'e'; 2-'f')
        # N - ranking number of energy levels of (p,j,c) blocks
        # n - total rank of (p,j,c) blocks
        #
        # Calculation limitations:
        # pmax = 40
        # jmax = 300
        # Ecut = 44600 cm-1
        # ---------------
        p   c   j   N   n   E   Evib    Erot    jref
        0   1   0   1   1   0.000   0.000   0.000   0
        0   1   2   1   1   2.341   0.000   2.341   0
        0   1   4   1   1   7.804   0.000   7.804   0
        0   1   6   1   1   16.389  0.000   16.389  0
        0   1   8   1   1   28.095  0.000   28.095  0

    See an example in `test/files/co2_cdsd_hamiltonian_fragment.levels <https://github.com/radis/radis/blob/develop/radis/test/files/cdsd_hitemp_09_fragment.txt>`

    See Also
    --------

    :py:class:`~radis.levels.partfunc_cdsd.PartFuncCO2_CDSDtab`
    """

    def __init__(
        self,
        energy_levels,
        isotope,
        levelsfmt,  # ='cdsd-pc',
        use_cached=True,
        use_json=None,  # TODO: Deprecated, remove
        verbose=True,
        mode="full summation",
    ):

        # %% Init

        # Initialize PartitionFunctionCalculator for this electronic state
        ElecState = ElectronicState("CO2", isotope, "X", "1Σu+")
        super(PartFuncCO2_CDSDcalc, self).__init__(ElecState, mode=mode)

        # Check inputs ('return' is not mentionned in signature. it will just return
        # after cache name is given)
        assert use_cached in [True, False, "regen", "force", "return"]
        if isotope not in [1, 2]:
            raise ValueError(
                "CDSD Energies not defined for isotope: {0}".format(isotope)
            )
        if use_json is not None:
            warn(
                DeprecationWarning(
                    "use_json replaced with faster HDF5-based use_cached"
                )
            )
        # Get vibrational level definitions that match Energy Database (in particular
        # how Evib and Erot are calculated)
        # This is needed to be able to match the levels in the Line Database and
        # the levels in the Energy database
        if levelsfmt == "cdsd-p":
            viblvl_label = "p"
        elif levelsfmt == "cdsd-pc":
            viblvl_label = "pc"
        elif levelsfmt == "cdsd-pcN":
            viblvl_label = "pcN"
        elif levelsfmt == "cdsd-hamil":
            viblvl_label = "pcJN"
        elif levelsfmt is None:
            # dont label the levels. Wont be able to use the EnergyDatabase to fetch
            # vibrational energies for lines, however it can still be used to
            # calculate Partition functions independently from a Spectrum calculation
            viblvl_label = None
        else:
            raise ValueError(
                "Unknown Energy database format: levelsfmt = `{0}`".format(levelsfmt)
                + ". Use one of: `cdsd-p`, `cdsd-pc`, `cdsd-pcN`,`cdsd-hamil`"
            )

        # Store defaults
        self.verbose = verbose
        self.use_cached = use_cached
        self.levelsfmt = levelsfmt
        self.viblvl_label = viblvl_label
        self.last_modification = time.ctime(
            getmtime(getTestFile(r"co2_cdsd_hamiltonian_fragment.levels"))
        )
        if verbose >= 2:
            print("Last modification time: {0}".format(self.last_modification))

        # Get variables to store in metadata  (after default values have been set)
        molecule = "CO2"  # will be stored in cache file metadata
        last_modification = time.ctime(
            getmtime(getTestFile(r"co2_cdsd_hamiltonian_fragment.levels"))
        )

        _discard = [
            "self",
            "energy_levels",
            "verbose",
            "ElecState",
            "electronic_state",
            "use_json",
            "use_cached",
        ]
        # (dev) locals() automatically stores all variables: levelsfmt, viblvl_label, etc.
        metadata = filter_metadata(locals(), discard_variables=_discard)

        # %% Get levels
        # Function of use_cached value:
        # ... if True, use (and generate if doesnt exist) cache file.
        # ... if 'regen', regenerate cache file. If 'force', raise an error
        # ... if file doesnt exist.
        # If file is deprecated, regenerate it unless 'force' was used

        # Load cache file if exists

        cachefile = energy_levels + ".h5"
        self.cachefile = cachefile

        # If return, return after cachefile generated (used for tests)
        if use_cached == "return":
            return

        df = load_h5_cache_file(
            cachefile,
            use_cached,
            valid_if_metadata_is=metadata,
            relevant_if_metadata_above={},
            relevant_if_metadata_below={},
            current_version=radis.__version__,
            last_compatible_version=radis.config["OLDEST_COMPATIBLE_VERSION"],
            verbose=verbose,
        )

        if df is None:  # Read normal file
            df = pd.read_csv(energy_levels, comment="#", delim_whitespace=True)
            df = self._add_degeneracies(df)
            df = self._add_levels(df)

        self.df = df  # Store

        if use_cached and not exists(cachefile):
            save_to_hdf(
                self.df,
                cachefile,
                metadata=metadata,
                version=radis.__version__,
                key="df",
                overwrite=True,
                verbose=verbose,
            )

    def _add_degeneracies(self, df):
        """Calculate and store degeneracies in database df.

        Parameters
        ----------

        df: pandas Dataframe
            energy database

        Notes
        -----

        .. warning::

            we use the same energies as CO2 626 but with gi = 2 for CO2 isotope 2 (636).
            It is done in the __init__ method
        """
        # Rotational degeneracy
        gj = 2 * df.j + 1
        # ... state dependant (forbidden rotational level are not in the database):
        gs = 1
        # ... state independant:
        gi = self.gi()
        grot = gj * gs * gi

        # Vibrational degeneracy
        gvib = 1  # CDSD is rovibrational complete

        # Store
        df["gvib"] = gvib
        df["gj"] = gj
        df["grot"] = grot

        return df

    def _add_levels(self, df):

        viblvl_label = self.viblvl_label

        if viblvl_label == "p":
            df["viblvl"] = vib_lvl_name_cdsd_p(
                df.p,
            )
        elif viblvl_label == "pc":
            df["viblvl"] = vib_lvl_name_cdsd_pc(df.p, df.c)
        elif viblvl_label == "pcN":
            df["viblvl"] = vib_lvl_name_cdsd_pcN(df.p, df.c, df.N)
        elif viblvl_label == "pcJN":
            df["viblvl"] = vib_lvl_name_cdsd_pcJN(df.p, df.c, df.j, df.N)
        elif viblvl_label is None:
            # dont label the levels. Wont be able to use the EnergyDatabase to fetch
            # vibrational energies for lines, however it can still be used to
            # calculate Partition functions independently from a Spectrum calculation
            pass
        else:
            raise ValueError("Unexpected viblvl_label value: {0}".format(viblvl_label))

        return df

    def gs(self):
        from radis.db.degeneracies import gs

        I = self.isotope
        return gs(2, I)

    def gi(self):
        from radis.db.degeneracies import gi

        I = self.isotope
        return gi(2, I)


# %% Test
if __name__ == "__main__":

    from radis.test.levels.test_partfunc import _run_testcases

    print("Testing parfunc: {0}".format(_run_testcases()))
