# -*- coding: utf-8 -*-
"""Partition function calculators and tabulators.

Notes
-----

Partition function calculators and tabulators

Calculators all derive from the same RovibPartitionFunction object,
and require a list of energies

Tabulators are more specific, and require a list of precalculated
partition functions at different temperature. PartFuncHAPI
uses the HITRAN hapi.py [1]_ module to interpolate Q for HITRAN species


Routine Listing
---------------

Partition functions:

- :class:`~radis.levels.partfunc.PartFuncHAPI`
- :class:`~radis.levels.partfunc.PartFunc_Dunham`

Which inherit from:

- :class:`~radis.levels.partfunc.RovibParFuncCalculator`
- :class:`~radis.levels.partfunc.RovibParFuncTabulator`

Which inherit from:

- :class:`~radis.levels.partfunc.RovibPartitionFunction`


See Also
--------

- :class:`~radis.levels.partfunc_cdsd.PartFuncCO2_CDSDcalc`
- :class:`~radis.levels.partfunc_cdsd.PartFuncCO2_CDSDtab`

References
----------

.. [1] `HAPI: The HITRAN Application Programming Interface <http://hitran.org/hapi>`_



-------------------------------------------------------------------------------
"""

# TODO: vectorize partition function caclulations for different temperatures. Would need
# stuff like E = df.E.values.reshape((1,-1)), etc.

# TODO: store molecule_data.json in the H5 file metadata. If not done already.


import sys
from os.path import exists
from warnings import warn

import numpy as np
import pandas as pd
from numpy import exp

import radis
from radis import OLDEST_COMPATIBLE_VERSION
from radis.db.classes import (
    HITRAN_CLASS1,
    HITRAN_CLASS2,
    HITRAN_CLASS3,
    HITRAN_CLASS5,
    HITRAN_CLASS6,
    get_molecule_identifier,
)
from radis.io.cache_files import load_h5_cache_file, save_to_hdf
from radis.lbl.labels import vib_lvl_name_hitran_class1, vib_lvl_name_hitran_class5
from radis.misc.basics import all_in
from radis.misc.debug import printdbg
from radis.misc.progress_bar import ProgressBar
from radis.phys.constants import hc_k  # ~ 1.44 cm.K


class RovibPartitionFunction(object):
    """General class from which all partition function calculators derive.

    Parameters
    ----------
    electronic_state: :class:`~radis.db.classes.ElectronicState`
        an :class:`~radis.db.classes.ElectronicState` object, which is
        defined in RADIS molecule database and contains spectroscopic data

    Notes
    -----
    Implementation:

    one partition function generator (RovibPartitionFunction) is generated
    per specie per isotope

    RovibPartitionFunction may differ in the way they build / fetch
    their list of states and the associated energies, but
    the .at(),  .at_noneq()  calls to calculate the partition
    function should be shared among all derived classes

    See Also
    --------
    :class:`~radis.levels.partfunc.RovibParFuncTabulator`,
    :class:`~radis.levels.partfunc.RovibParFuncCalculator`
    """

    def __init__(self, electronic_state):

        #        # Check inputs
        #        if type(molecule) is int:
        #            molecule = get_molecule(molecule)
        #
        #        # Get molecule from radis.db
        #        mol_list = import_from_module('radis.db.molecules', molecule)
        #        try:
        #            mol = mol_list[isotope]
        #        except KeyError:
        #            raise KeyError('Isotope {0} not defined for molecule {1}'.format(int(isotope), molecule))

        ElecState = electronic_state

        try:
            ElecState.Erovib
        except AttributeError:
            raise AttributeError(
                "{0} has no energy function defined in RADIS".format(
                    ElecState.get_fullname()
                )
            )

        # Store
        self.ElecState = ElecState  # electronic state object
        self.molecule = ElecState.name  # molecule name
        self.isotope = ElecState.iso

        # Line database
        # pandas Dataframe that holds all the lines
        self.df = pd.DataFrame({})
        # updated on inherited classes initialization


# %% Subclasses :
# - Tabulator
# - Calculator


class RovibParFuncTabulator(RovibPartitionFunction):
    def at(self, T, **kwargs):
        """Get partition function at temperature T under equilibrium
        conditions, from tabulated data.

        Parameters
        ----------

        T: float
            equilibrium temperature

        Returns
        -------

        Q: float
            partition function interpolated  at temperature T

        Examples
        --------

        See :ref:`online examples <label_examples_partition_functions>`

        See Also
        --------

        :py:class:`~radis.levels.partfunc.PartFuncHAPI`
        """

        # For compatibility with the syntax of RovibParFuncCalculator.at()
        update_populations = kwargs.pop("update_populations", False)
        if kwargs != {}:
            raise TypeError(
                "RovibParFuncTabulator.at() got an unexpected keyword: {0}".format(
                    kwargs.keys()
                )
            )
        if update_populations != False:
            raise ValueError(
                "Cannot use update_populations with tabulated partition functions. "
                + "Use a direct summation instead. See RovibParFuncCalculator"
            )

        # defined individually for each class Variants (one per database)
        return self._at(T)

    def at_noneq(self, *args, **kwargs):
        raise ValueError(
            "Cannot calculate non equilibrium partition "
            + "functions from (equilibrium) tabulated values. Use a "
            + "Partition Function Calculator"
        )


class RovibParFuncCalculator(RovibPartitionFunction):
    """
    Parameters
    ----------

    electronic_state: :class:`~radis.db.classes.ElectronicState`
        an :class:`~radis.db.classes.ElectronicState` object, which is
        defined in RADIS molecule database and contains spectroscopic data.
        Energies are calculated with the :py:meth:`~radis.db.classes.ElectronicState.Erovib`
        method.

    See Also
    --------

    :py:class:`~radis.levels.partfunc.PartFunc_Dunham`,
    :py:meth:`~radis.db.classes.ElectronicState.Erovib`
    """

    def __init__(self, electronic_state):

        super(RovibParFuncCalculator, self).__init__(electronic_state=electronic_state)

    def at(self, T, update_populations=False):
        """Get partition function at temperature T under equilibrium
        conditions.

        Parameters
        ----------

        T: float
            equilibrium temperature

        Other Parameters
        ----------------

        update_populations: boolean
            if ``True``, store calculated populations in energy level list
            Default ``False``

        Returns
        -------

        Q: float
            partition function calculated at temperature T

        Examples
        --------

        See :ref:`online examples <label_examples_partition_functions>`

        See Also
        --------

        :meth:`~radis.levels.partfunc.RovibPartitionFunction.at_noneq`,
        :meth:`~radis.levels.partfunc.RovibPartitionFunction.at_noneq_3Tvib`
        """
        if __debug__:
            printdbg(
                "called RovibPartitionFunction.at(T={0}K, ".format(T)
                + "update_populations={0})".format(update_populations)
            )

        # Check inputs, initialize
        assert isinstance(update_populations, bool)
        df = self.df
        if "g" in df.columns:
            g = df.g
        elif all_in(["gvib", "grot"], df.columns):
            g = df.gvib * df.grot
        else:
            raise ValueError(
                "either g, or gvib+grot must be defined to "
                + "calculate total degeneracy. Got: {0}".format(list(df.keys()))
            )

        # Calculate

        nQ = g * exp(-hc_k * df.E / T)
        Q = nQ.sum()

        # Update energy level table with populations (doesnt
        # cost much and can be used to plot populations afterwards)
        # ... 'n'
        if update_populations:
            df["n"] = nQ / Q

        return Q

    def at_noneq(
        self,
        Tvib,
        Trot,
        overpopulation=None,
        vib_distribution="boltzmann",
        rot_distribution="boltzmann",
        returnQvibQrot=False,
        update_populations=False,
    ):
        """Calculate Partition Function under non equilibrium (Tvib, Trot),
        with boltzmann/treanor distributions and overpopulations as specified
        by the user.

        Parameters
        ----------

        Tvib, Trot: float

        overpopulation: dict, or ``None``
            dict of overpopulated levels: ``{'level':over_factor}``

        vib_distribution: ``'boltzmann'``, ``'treanor'``
            distribution of vibrational levels

        rot_distribution: ``'boltzmann'``
            distribution of rotational levels

        returnQvibQrot: boolean
            cf output

        Other Parameters
        ----------------

        update_populations: boolean
            if ``True``, store calculated populations in energy level list.
            Default ``False``

        Returns
        -------

        Q: float
            partition function calculated at non eq temperatures

        if returnQvibQrot:

        Q, Qvib, dfQrot: float, float, pandas table
            total partition function, vibrational partition function,
            and table of rotational partition functions for each vibrational
            state (note that all Qrot are not necessarily the same
            for all vibrational levels)

        Examples
        --------

        See :ref:`online examples <label_examples_partition_functions>`

        See Also
        --------

        :meth:`~radis.levels.partfunc.RovibPartitionFunction.at`,
        :meth:`~radis.levels.partfunc.RovibPartitionFunction.at_noneq_3Tvib`
        """
        if __debug__:
            printdbg(
                "called RovibPartitionFunction.atnoneq"
                + "(Tvib={0}K, Trot={1}K, ... )".format(Tvib, Trot)
                + "update_populations={0})".format(update_populations)
            )
        #                               'overpopulation={0}, vib_distribution={1}'.format(overpopulation, vib_distribution)+\
        #                               'rot_distribution={0}'.format(rot_distribution)

        # Check inputs, initialize
        if overpopulation is None:
            overpopulation = {}
        else:
            if not returnQvibQrot:
                raise ValueError(
                    "When using overpopulation, partition function "
                    + "must be calculated with returnQvibQrot=True"
                )
        assert vib_distribution in ["boltzmann", "treanor"]
        assert rot_distribution in ["boltzmann"]
        if vib_distribution == "boltzmann":
            if not "Evib" in list(self.df.keys()):
                raise ValueError(
                    "Evib must be defined to calculate non-equilibrium "
                    + "partition functions"
                )
        elif vib_distribution == "treanor":
            if not all_in(["Evib_h", "Evib_a"], list(self.df.keys())):
                raise ValueError(
                    "Evib_h and Evib_a must be defined to calculate non-equilibrium "
                    + "partition functions with treanor distribution"
                )
        if rot_distribution == "boltzmann":
            if not "Erot" in list(self.df.keys()):
                raise ValueError(
                    "Evib and Erot must be defined to calculate non-equilibrium "
                    + "partition functions"
                )

        # Get variables
        df = self.df
        gvib = df.gvib  # self.gvib(M, I)
        grot = df.grot
        # Calculate

        # ... mode: Trot, Tvib, + overpopulation
        if returnQvibQrot:

            if not "viblvl" in self.df:
                raise KeyError(
                    "To return Qrot and Qvib vibrational levels must be "
                    + "identified in the database. Lookup function "
                    + "add_bands in radis.lbl.bands"
                )

            # ... Vibrational populations
            if vib_distribution == "boltzmann":
                df["nvibQvib"] = gvib * exp(-df.Evib * hc_k / Tvib)
            elif vib_distribution == "treanor":
                df["nvibQvib"] = gvib * exp(
                    -hc_k * (df.Evib_h / Tvib + df.Evib_a / Trot)
                )
            else:
                raise NotImplementedError
            # ... Rotational populations
            if rot_distribution == "boltzmann":
                df["nrotQrot"] = grot * exp(-df.Erot * hc_k / Trot)
            else:
                raise NotImplementedError

            # If overpopulation, first check all levels exist
            levels = df.viblvl.unique()
            for viblvl in overpopulation.keys():
                if not viblvl in levels:
                    raise ValueError(
                        "Level {0} not in energy levels database".format(viblvl)
                    )
                    # could be a simple warning too
            # Add overpopulations (so they are taken into account in the partition function)
            for viblvl, ov in overpopulation.items():
                if ov != 1:
                    df.loc[df.viblvl == viblvl, "nvibQvib"] *= ov
                    # TODO: add warning if empty? I dont know how to do it without
                    # an extra lookup though.

            # Calculate sum of levels
            nQ = df.nvibQvib * df.nrotQrot
            Q = nQ.sum()

            # Group by vibrational level, and get level-dependant
            # quantities such as vib degeneracy, Evib, etc.
            dgb = df.groupby(["viblvl"], as_index=True)
            # perf: 345 ms ± 5.88 ms on test case
            viblvl_Qrot = dgb["nrotQrot"].sum()
            # (all the same in one vib group)
            viblvl_nvibQvib = dgb["nvibQvib"].first()
            Qvib = viblvl_nvibQvib.sum()  # (here we only sum once per vib group)
            # (all the same in one vib group)
            viblvl_Evib = dgb["Evib"].first()
            # (all the same in one vib group)
            viblvl_gvib = dgb["gvib"].first()

            # Energies, degeneracies, populations for each vibrational band
            dfQrot = pd.DataFrame({"Qrot": viblvl_Qrot})
            dfQrot["nvib"] = viblvl_nvibQvib / Qvib
            dfQrot["Evib"] = viblvl_Evib
            dfQrot["gvib"] = viblvl_gvib

            # Update energy level table with populations (doesnt
            # cost much and can be used to plot populations afterwards)
            # Adds: 'nvib', 'n', 'nrot', 'Qrot'
            if update_populations:
                df["nvib"] = df.nvibQvib / Qvib
                df["n"] = nQ / Q

                # get rotational populations
                # ... reindexing dfQrot to get a direct access by viblvl
                dfQrot_dict = dict(list(zip(dfQrot.index, dfQrot.Qrot)))

                # ... Add Qrot
                df_viblvl = df.set_index(["viblvl"], inplace=False)
                df["Qrot"] = df_viblvl.index.map(dfQrot_dict.get).values

                # ... infer nrot
                df["nrot"] = df.nrotQrot / df.Qrot

            # Check that partition functions are valid
            if __debug__:  # discarded if running with python -O
                # the calculation below is valid without Born-Oppenheimer, however
                # it assumes Boltzmann distribution + no overpopulation
                #                Qrovib = ((dfQrot.gvib*exp(-dfQrot.Evib*hc_k/Tvib))*dfQrot.Qrot).sum()
                # general case:
                Qrovib = (dfQrot.nvib * Qvib * dfQrot.Qrot).sum()
                if not np.isclose(Q, Qrovib):
                    raise ValueError(
                        "Rovibrational partition function ({0:.2f}) doesnt ".format(Q)
                        + "match value recomputed from vibrational and rotational "
                        + "partition functions ({0:.2f}). Check how Evib and Erot ".format(
                            Qrovib
                        )
                        + "are defined in your Energy Database"
                    )

            # Clean
            del df["nrotQrot"]
            del df["nvibQvib"]

            return Q, Qvib, dfQrot

        # ... mode: Trot, Tvib, no overpopulation, no vib band details
        else:  # slightly faster, but doesnt return nvib nor Qvib

            g = gvib * grot

            if vib_distribution == "boltzmann" and rot_distribution == "boltzmann":
                nQ = g * exp(-hc_k * (df.Evib / Tvib + df.Erot / Trot))
            elif vib_distribution == "treanor" and rot_distribution == "boltzmann":
                nQ = g * exp(-hc_k * (df.Evib_h / Tvib + (df.Evib_a + df.Erot) / Trot))
            else:
                raise NotImplementedError

            Q = nQ.sum()

            # Update energy level table with populations (doesnt
            # cost much and can be used to plot populations afterwards)
            # ... add: 'n'
            if update_populations:
                df["n"] = nQ / Q

            return Q

    def at_noneq_3Tvib(
        self,
        Tvib,
        Trot,
        overpopulation=None,
        vib_distribution="boltzmann",
        rot_distribution="boltzmann",
        returnQvibQrot=False,
        update_populations=False,
    ):
        """Calculate Partition Function under non equilibrium ((Tvib1, Tvib2,
        Tvib3), Trot), with boltzmann/treanor distributions and overpopulations
        as specified by the user.

        Dedicated function for 3 Tvib mode

        Parameters
        ----------

        Tvib, Trot: float

        overpopulation: dict, or ``None``
            dict of overpopulated levels: ``{'level':over_factor}``

        vib_distribution: ``'boltzmann'``, ``'treanor'``

        rot_distribution: ``'boltzmann'``

        returnQvibQrot: boolean
            cf output

        update_populations: boolean
            if ``True``, saves populations for calculated temperature in PartitionFunction
            dataframe. Default ``False``

        Returns
        -------

        Q: float
            partition function calculated at non eq temperatures

        if returnQvibQrot:

        Q, Qvib, dfQrot: float, float, pandas table
            total partition function, vibrational partition function,
            and table of rotational partition functions for each vibrational
            state (note that all Qrot are not necessarily the same
            for all vibrational levels)

        See Also
        --------

        :meth:`~radis.levels.partfunc.RovibPartitionFunction.at`,
        :meth:`~radis.levels.partfunc.RovibPartitionFunction.at_noneq`
        """
        if __debug__:
            printdbg(
                "called RovibPartitionFunction.at_noneq_3Tvib"
                + "(Tvib={0}K, Trot={1}K)".format(Tvib, Trot)
                + "update_populations={0})".format(update_populations)
            )
        #                               'overpopulation={0}, vib_distribution={1}'.format(overpopulation, vib_distribution)+\
        #                               'rot_distribution={0})'.format(rot_distribution)

        # Check inputs
        if overpopulation is None:
            overpopulation = {}
        #        else:
        #            if not returnQvibQrot:
        #                raise ValueError('When using overpopulation partition function '+\
        #                                 'must be calculated with returnQvibQrot')
        assert vib_distribution in ["boltzmann", "treanor"]
        assert rot_distribution in ["boltzmann"]
        if vib_distribution == "boltzmann":
            if not all_in(["Evib1", "Evib2", "Evib3"], list(self.df.keys())):
                raise ValueError(
                    "Evib1, Evib2, Evib3 must be defined to calculate non-equilibrium "
                    + "partition functions with Tvib1, Tvib2, Tvib3. Got {0}".format(
                        list(self.df.keys())
                    )
                )
        elif vib_distribution == "treanor":
            if not all_in(
                ["Evib1_h", "Evib1_a", "Evib2_h", "Evib2_a", "Evib3_h", "Evib3_a"],
                list(self.df.keys()),
            ):
                raise ValueError(
                    "Evib1_h, Evib1_a, Evib2_h, Evib2_a, Evib3_h, Evib3_a "
                    + "must be defined to calculate non-equilibrium "
                    + "partition functions with treanor distribution "
                    + "and Tvib1, Tvib2, Tvib3"
                )
        if rot_distribution == "boltzmann":
            if not "Erot" in list(self.df.keys()):
                raise ValueError(
                    "Evib and Erot must be defined to calculate non-equilibrium "
                    + "partition functions"
                )

        # Get variables
        Tvib1, Tvib2, Tvib3 = Tvib
        df = self.df
        gvib = df.gvib
        grot = df.grot

        # Calculate
        if overpopulation:
            raise NotImplementedError("overpopulation not implemented for 3 Tvib mode")

        # ... mode: Trot, Tvib, + overpopulation
        if returnQvibQrot:
            raise NotImplementedError("returnQvibQrot not implemented for 3 Tvib mode")

        # ... mode: Trot, Tvib, no overpopulation, no vib band details
        else:  # much faster, but doesnt return nvib nor Qvib

            g = gvib * grot

            if vib_distribution == "boltzmann" and rot_distribution == "boltzmann":
                nQ = g * exp(
                    -hc_k
                    * (
                        df.Evib1 / Tvib1
                        + df.Evib2 / Tvib2
                        + df.Evib3 / Tvib3
                        + df.Erot / Trot
                    )
                )
            elif vib_distribution == "treanor" and rot_distribution == "boltzmann":
                nQ = g * exp(
                    -hc_k
                    * (
                        df.Evib1_h / Tvib1
                        + df.Evib2_h / Tvib2
                        + df.Evib3_h / Tvib3
                        + (df.Evib1_a + df.Evib2_a + df.Evib3_a + df.Erot) / Trot
                    )
                )
            else:
                raise NotImplementedError

            Q = nQ.sum()

            # Update energy level table with populations
            # ... add: 'n'
            if update_populations:
                df["n"] = nQ / Q

            return Q

    def reset_populations(self):
        """Discard computed populations of all energy levels.

        To call on every RovibrationalPartitionFunction object before
        each new spectrum calculation
        """

        for k in ["nvib", "n", "nrot"]:
            if k in self.df:
                del self.df[k]

    # %% Methods to get populations of all states

    def _get_vib_populations(self):
        """Return vibrational populations for all levels featured in given line
        set."""

        df = self.df

        return df.drop_duplicates("viblvl")

    def _get_rovib_populations(self):
        """Return rovibrational populations for all levels featured in the
        energy levels list df.

        Notes
        -----

        assumes a complete rovibrational assigmnent but no hyperfine assignment
        (i.e: all energy levels are returned!). If hyperfine assigmnent is given,
        this method should be modified to return only the ``roviblvl`` unique
        keys
        """

        df = self.df

        return df


# %% Variants of Tabulated partition functions (interpolate)


class PartFuncHAPI(RovibParFuncTabulator):
    """Return partition function using interpolation of tabulated values.

    Parameters
    ----------

    M: int
        molecule id

    I: int
        isotope identifier

    path: str
        path to ``hapi.py``. If None, the hapi package from PyPI is used.

    Examples
    --------

    ::

        from radis.levels.partfunc import PartFuncHAPI
        from radis.db.classes import get_molecule_identifier

        M = get_molecule_identifier('N2O')
        iso=1

        Q = PartFuncHAPI(M, iso)
        print(Q.at(T=1500))

    See :ref:`online examples <label_examples_partition_functions>` for more.

    References
    ----------

    Partition function are calculated with HAPI [1]_ (Hitran Python Interface) using
    partitionSum(M,I,T)

    .. [1] `HAPI: The HITRAN Application Programming Interface <http://hitran.org/hapi>`_
    """

    def __init__(self, M, I, path=None, verbose=True):

        self.verbose = verbose

        if path is not None:
            partitionSum = self.import_from_file(path)
        else:
            # Use RADIS embedded
            from hapi import HAPI_VERSION, partitionSum

            if self.verbose >= 2:
                print("HAPI version: %s" % HAPI_VERSION)

        # Check inputs
        if isinstance(M, str):
            M = get_molecule_identifier(M)
        if type(M) is not int:
            raise TypeError("Molecule id must be int: got {0} ({1})".format(M, type(M)))
        if type(I) is not int:
            raise TypeError(
                "Isotope number must be int: got {0} ({1})".format(I, type(I))
            )

        self.partitionSum = partitionSum
        self.M = M
        self.I = I

    def import_from_file(self, path):
        """Import hapi.py from a given file (in case user wants to specify a
        different HAPI version than the one embedded in RADIS)"""
        if sys.version == 2:
            import imp

            hapi = imp.load_source("hapi", path)
        else:  # Python 3
            import importlib.util

            spec = importlib.util.spec_from_file_location("hapi", path)
            hapi = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(hapi)
        return hapi.partitionSum

    def _at(self, T):
        """Get partition function of species M, isotope I at temperature T.

        Called by :meth:`radis.levels.partfunc.RovibParFuncTabulator.at`
        """
        return self.partitionSum(self.M, self.I, T)


# %% Calculated partition functions (either from energy levels, or ab initio)


def _get_cachefile_name(ElecState):
    """Get name of cache file for calculated rovibrational energies.

    Basically store it alongside the jsonfile of the ElecState with:

        [jsonfile]_[molecule_name]_[isotope]_[electronic state]_levels.h5
    """
    molecule = ElecState.name
    isotope = ElecState.iso
    state = ElecState.state
    jsonfile = ElecState.jsonfile
    filename = "{0}_{1}_iso{2}_{3}_levels.h5".format(
        jsonfile.replace(".json", ""), molecule, isotope, state
    )
    return filename


class PartFunc_Dunham(RovibParFuncCalculator):
    """Calculate partition functions from spectroscopic constants, if
    molecule data is available in RADIS. If using the
    :ref:`default spectroscopic constants <label_db_spectroscopic_constants>`,
    make sure you know what reference data is being used in RADIS!
    Energies are calculated with the :py:meth:`~radis.db.classes.ElectronicState.Erovib`
    method of :py:data:`~radis.db.molecules.Molecules`

    Parameters
    ----------

    electronic_state: :class:`~radis.db.classes.ElectronicState`
        an :class:`~radis.db.classes.ElectronicState` object, which is
        defined in RADIS molecule database and contains spectroscopic data

    vmax: int, or ``None``
        maximum vibrational quantum to calculate with Dunham expansion.
        If None, the molecule one is taken.
        If None still, all levels are calculated up to the molecule dissociation
        energy

    vmax_morse: int, or ``None``
        maximum vibrational quantum to calculate with Morse potential.
        If None, the molecule one is taken. Use ``0`` or ``-1`` not to calculate
        with Morse potential

    Jmax: int, or ``None``
        maximum rotational quantum. If None, the molecule one is taken.
        If None, all levels are calculated up to the molecule dissociation
        energy

    use_cached: ``True``, ``False``, or ``'regen'``, ``'force'``
        if ``True``, use (and generate if doesnt exist) a ``.h5`` file.
        If ``'regen'``, regenerate cache file. If ``'force'``, raise an error
        if file doesnt exist. Default ``True``

    Other Parameters
    ----------------

    calc_Evib_per_mode: boolean
        if ``True``, calculate energies of each vibration mode (so far only
        implemented for CO2 with Evib1, Evib2, Evib3 but shall be generalized
        to all molecules)

    calc_Evib_harmonic_anharmonic: boolean
        if ``True``, calculate and store separately harmonic and anharmonic parts
        of the vibrational energy. This is needed to calculate Treanor distributions
        ( ~ Evib_harmonic / Tvib  - Evib_anharmonic / Trot )

    group_energy_modes_in_2T_model: dict
        (experimental in neq 0.9.21) for polyatomic molecules (> 1 vibration mode),
        how to regroup energy modes when working with 2T models. For instance,
        for CO2, (Evib, Erot) could as well as defined with::

            ['Evib1', 'Evib2', 'Evib3'],['Erot']
            or
            ['Evib3'],['Evib1', 'Evib2', 'Erot']

        depending on which levels are supposed to interact the most

    Notes
    -----

    Validity:

    So far, RADIS energy levels are only calculated from Dunham's expansion.
    Above a certain vibrational level a Morse Potential may be used. See how
    the molecule is defined in :mod:`~radis.db.molecules`

    See Also
    --------

    :mod:`~radis.db.molecules`,
    :py:meth:`~radis.db.classes.ElectronicState.Erovib`,
    :py:data:`~radis.db.molecules.Molecules`

    """

    def __init__(
        self,
        electronic_state,
        vmax=None,
        vmax_morse=None,
        Jmax=None,
        spectroscopic_constants="default",
        use_cached=True,
        verbose=True,
        calc_Evib_per_mode=True,
        calc_Evib_harmonic_anharmonic=True,
        group_energy_modes_in_2T_model={"CO2": (["Evib1", "Evib2", "Evib3"], ["Erot"])},
    ):  # , ZPE=None):
        # TODO: find a way to initialize calc_Evib_per_mode or calc_Evib_harmonic_anharmonic from
        # the SpectrumFactory on Spectrum generation time...
        # Maybe recompute the cache file if needed?

        # TODO @dev: refactor : move group_energy_modes_in_2T_model in radis/config.json ?

        # WIP @erwan: spectroscopic_constants: NOT USED for the moment
        # parameter to change the constants
        # used to calculate the electronic_state energies. An alternative strategy
        # be to adjust the constants in electronic_state directly.

        # %% Init
        super(PartFunc_Dunham, self).__init__(electronic_state=electronic_state)

        # Check inputs ('return' is not mentionned in signature. it will just return
        # after cache name is given)
        assert use_cached in [True, False, "regen", "force", "return"]

        ElecState = self.ElecState  # molecule
        molecule = ElecState.name  # will be stored in cache file metadata
        isotope = ElecState.iso  # will be stored in cache file metadata
        if vmax is None:
            vmax = ElecState.vmax
        if vmax_morse is None:
            vmax_morse = ElecState.vmax_morse
        if Jmax is None:
            Jmax = ElecState.Jmax

        # How to lump energy modes
        if molecule in group_energy_modes_in_2T_model:
            group_energy_modes = group_energy_modes_in_2T_model[molecule]

        # Store
        self.vmax = vmax
        self.vmax_morse = vmax_morse
        self.Jmax = Jmax
        self.verbose = verbose
        self.use_cached = use_cached
        self.group_energy_modes_in_2T_model = group_energy_modes_in_2T_model

        # Get variables to store in metadata  (after default values have been set)
        # ... this ensures that cache files generated with different metadata
        # ... will not be used (in particular: different rovib_constants)
        metadata = {
            "molecule": molecule,
            "isotope": isotope,
            "rovib_constants": ElecState.rovib_constants,
            "vmax": vmax,
            "vmax_morse": vmax_morse,
            "Jmax": Jmax,
            "group_energy_modes_in_2T_model": group_energy_modes_in_2T_model,
            "calc_Evib_harmonic_anharmonic": calc_Evib_harmonic_anharmonic,
            "calc_Evib_per_mode": calc_Evib_per_mode,
        }
        if molecule in group_energy_modes_in_2T_model:
            metadata.update({"group_energy_modes": group_energy_modes})

        # get cache file path

        # Function of use_cached value:
        # ... if True, use (and generate if doesnt exist) cache file.
        # ... if 'regen', regenerate cache file. If 'force', raise an error
        # ... if file doesnt exist.
        # If file is deprecated, regenerate it unless 'force' was used

        cachefile = _get_cachefile_name(ElecState)
        self.cachefile = cachefile

        # If return, return after cachefile generated (used for tests)
        if use_cached == "return":
            return

        self.df = load_h5_cache_file(
            cachefile,
            use_cached,
            valid_if_metadata_is=metadata,
            relevant_if_metadata_above={},
            relevant_if_metadata_below={},
            current_version=radis.__version__,
            last_compatible_version=OLDEST_COMPATIBLE_VERSION,
            verbose=verbose,
        )

        # Get levels
        if molecule in HITRAN_CLASS1 + HITRAN_CLASS2 + HITRAN_CLASS3:
            nb_vib_modes = 1
        elif molecule in HITRAN_CLASS5 + HITRAN_CLASS6:
            nb_vib_modes = 3
        else:
            raise NotImplementedError
        self.nb_vib_modes = nb_vib_modes

        # Build energy levels if needed
        if self.df is None:
            # Build energy levels
            if verbose:
                print(
                    "Calculating energy levels with Dunham expansion for {0}".format(
                        ElecState.get_fullname()
                    )
                )
                if not use_cached:
                    print(
                        "Tip: set ``use_cached=True`` next time not to recompute levels every time"
                    )
            if molecule in HITRAN_CLASS1:
                self.build_energy_levels_class1()
            elif molecule in HITRAN_CLASS5:  # CO2
                self.build_energy_levels_class5(
                    calc_Evib_per_mode=calc_Evib_per_mode,
                    calc_Evib_harmonic_anharmonic=calc_Evib_harmonic_anharmonic,
                    group_energy_modes_in_2T_model=group_energy_modes,
                )
            else:
                raise NotImplementedError

        # save files if use_cached
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

        # Add extra columns (note that this is not saved to disk)
        self._add_extra()

        return

    def _add_extra(self):

        df = self.df
        molecule = self.molecule

        if molecule in ["CO2"]:
            # Add harmonic / anharmonic energies
            if (
                all_in(
                    ["Evib1_a", "Evib1_h", "Evib2_a", "Evib2_h", "Evib3_a", "Evib3_h"],
                    df,
                )
                and not "Evib_h" in df
                and not "Evib_a" in df
            ):

                group_energy_modes = self.group_energy_modes_in_2T_model[molecule]

                # Add total harmonic / Anharmonic energies as the sum of each mode
                # (used in 2 T temperature models)

                if group_energy_modes == (["Evib1", "Evib2", "Evib3"], ["Erot"]):
                    df["Evib_h"] = df.Evib1_h + df.Evib2_h + df.Evib3_h
                    df["Evib_a"] = df.Evib1_a + df.Evib2_a + df.Evib3_a

                elif group_energy_modes == (["Evib3"], ["Evib1", "Evib2", "Erot"]):
                    if not all_in(["Evib1", "Evib2", "Evib3"], df):
                        raise KeyError(
                            "You need Evib1, Evib2, Evib3 calculated separately. "
                            + "Use calc_Evib_per_mode=True"
                        )
                    df["Evib_h"] = df.Evib3_h
                    df["Evib_a"] = df.Evib3_a

                else:
                    raise NotImplementedError(
                        "group_energy_modes_in_2T_model: {0}".format(group_energy_modes)
                    )

    def build_energy_levels_class1(self):  # , ZPE=0):
        """in the case where only Ediss is given. Deal with vmax, Jmax later.

        Applies to molecules in :data:`~radis.io.hitran.HITRAN_CLASS1`

        Returns
        -------

        None:
            but the Pandas DataFrame `self.df` is updated with parameters:

            - ``g`` : degeneracy
            - ``E`` : energy level
            - ``Evib`` : vibrational energy
            - ``Erot`` : rotational energy
            - ``viblvl`` : vibrational level name
        """

        vib_lvl_name = vib_lvl_name_hitran_class1

        ElecState = self.ElecState  # molecule

        Ediss = ElecState.Ediss  # cm-1

        vmax = self.vmax  # max vibrational level to calculate with Dunham
        # max vibrational level to calculate with Morse potential
        vmax_morse = self.vmax_morse
        Jmax = self.Jmax  # max rotational level to calculate

        # First get max vib level  (to be calculated with Dunham expansion)
        if vmax is None:
            vmax = 0
            while ElecState.Erovib(vmax, J=0, remove_ZPE=True) < Ediss:
                vmax += 1
        if Jmax is None:
            Jmax = np.nan  # no limit
        Jmaxcalc = 0  # max calculated. just for info

        # Calculate lower rovibrational levels with Dunham expansion
        # ------------
        # format: (v, j, E, Evib)   # g, Erot calculated at the end (faster)
        levels = []
        Evib = ElecState.Erovib(v=0, J=0, remove_ZPE=True)  # no Zero-point-energy

        # ... Loop over vibrational levels:
        for v in range(0, vmax + 1):
            Evib = ElecState.Erovib(v, J=0, remove_ZPE=True)  # no Zero-point-energy
            if __debug__:
                printdbg(
                    "Calculating Evib for "
                    + "v={0}: {1:.2f}cm-1 (Dunham expansion)".format(v, Evib)
                )
            # ... Rotational loop
            E = Evib
            J = 0
            while 0 <= E < Ediss and not J > Jmax:  # (passes if Jmax is nan)
                # store rovib level
                levels.append([v, J, E, Evib])
                # calculate new one:
                J += 1
                E = ElecState.Erovib(v, J, remove_ZPE=True)  # no Zero-point-energy
            Jmaxcalc = max(Jmaxcalc, J)

        # If defined, calculate upper rovibrational levels with Morse potential
        # ----------
        if vmax_morse is not None:

            # vibrational energy of the last vibrational level for which Dunham expansion is valid
            Evib_last = ElecState.Erovib(vmax, 0)
            # difference of the 2 last vib level for which Dunham expansion is valid
            delta_E_last = Evib_last - ElecState.Erovib(vmax - 1, 0)

            v_inc = ElecState.get_Morse_inc()

            # ... Start loop on all Morse potential levels
            Evib = Evib_last
            for v in range(vmax + 1, vmax_morse + 1):
                delta_E = delta_E_last - (v + 1 - vmax) * v_inc
                Evib = Evib + delta_E
                if Evib > Ediss:
                    warn(
                        "Energy above dissociation threshold: {0}>{1}".format(
                            Evib, Ediss
                        )
                    )
                if __debug__:
                    printdbg(
                        "Calculating Evib for "
                        + "v={0}: {1:.2f}cm-1 (Morse Potential)".format(v, Evib)
                    )
                # ... Rotational loop
                J = 0
                E = Evib
                # (passes if Jmax is nan)
                while 0 <= E < Ediss and not J >= Jmax:
                    # store rovib level
                    levels.append([v, J, E, Evib])
                    # calculate new one:
                    J += 1
                    # no Zero-point-energy
                    Erot = ElecState.Erovib(0, J, remove_ZPE=True)
                    E = Evib + Erot
                Jmaxcalc = max(Jmaxcalc, J)

        df = pd.DataFrame(levels, columns=["v", "j", "E", "Evib"])

        # Store vibrational level name
        df["viblvl"] = vib_lvl_name(df.v)

        # Calculate degeneracies
        df["gj"] = 2 * df.j + 1
        df["gvib"] = 1  # energy base is assumed to be rovibrational complete
        gs = self.gs(ElecState)

        if isinstance(gs, tuple):
            raise NotImplementedError(
                "Different degeneracies (symmetric/antisymmetric)"
                + " not implemented for molecule {0} isotope {1}".format(
                    ElecState.id, ElecState.iso
                )
            )
            # (for devs: to implement that efficiently see CO2 (_class5))

        gi = self.gi(ElecState)
        df["grot"] = gs * gi * df.gj
        df["Erot"] = df.E - df.Evib

        if self.verbose:
            print("Database generated up to v={0}, J={1}".format(v, Jmaxcalc))

        self.df = df

    def build_energy_levels_class5(
        self,
        calc_Evib_per_mode=True,
        calc_Evib_harmonic_anharmonic=False,
        group_energy_modes_in_2T_model=(["Evib1", "Evib2", "Evib2"], ["Erot"]),
    ):
        """in the case where only Ediss is given. Deal with vmax, Jmax later.

        :data:`~radis.io.hitran.HITRAN_CLASS5` = ['CO2']
        # Linear triatomic with large Fermi resonance

        Parameters
        ----------

        calc_Evib_per_mode: boolean
            if ``True``, calculates Evib1, Evib2, Evib3

        Other Parameters
        ----------------

        calc_Evib_harmonic_anharmonic: boolean
            if ``True``, calculates Evib_h and Evib_a, harmonic and non harmonic
            components, to be used in Treanor distributions

        group_energy_modes_in_2T_model: dict
            (experimental in neq 0.9.21) for polyatomic molecules (> 1 vibration mode),
            how to regroup energy modes when working with 2T models. For instance,
            for CO2, (Evib, Erot) could as well as defined with::

                ['Evib1', 'Evib2', 'Evib3'],['Erot']

            or::

                ['Evib3'],['Evib1', 'Evib2', 'Erot']

            depending on which levels are supposed to interact the most

        Returns
        -------

        None:
            but the Pandas DataFrame `self.df` is updated with parameters:

            - ``g`` : degeneracy
            - ``E`` : energy level
            - ``Evib`` : vibrational energy
            - ``Erot`` : rotational energy
            - ``viblvl`` : vibrational level name
        """

        vib_lvl_name = vib_lvl_name_hitran_class5

        ElecState = self.ElecState
        Ediss = ElecState.Ediss  # cm-1

        vmax = self.vmax
        vmax_morse = self.vmax_morse
        Jmax = self.Jmax
        if vmax_morse is not None:
            raise NotImplementedError

        if vmax is None:
            vmax = ElecState.vmax
        if Jmax is None:
            Jmax = ElecState.Jmax

        #        if vmax is None: vmax = 100    # just to prevent infinite loop. Ediss should be limitant
        if Jmax is None:
            Jmax = 1000

        # First get max levels
        v1max = 0
        v2max = 0
        v3max = 0
        while ElecState.Erovib(v1max, 0, 0, 0, J=0, remove_ZPE=True) < Ediss:
            v1max += 1
        while ElecState.Erovib(0, v2max, v2max, 0, J=0, remove_ZPE=True) < Ediss:
            v2max += 1
        while ElecState.Erovib(0, 0, 0, v3max, J=0, remove_ZPE=True) < Ediss:
            v3max += 1
        if vmax is not None:
            v1max = min(vmax, v1max)
            v2max = min(vmax, v2max)
            v3max = min(vmax, v3max)

        # Then fill mixed modes levels

        # %%

        levels = []  # vibrational levels
        pb = ProgressBar(v1max + 1, active=True)
        Jmax_calc = 0
        for v1 in range(v1max + 1):
            pb.update(v1)
            for v2 in range(v2max + 1):
                for l2 in [v2]:
                    # Calculation with HITRAN spectroscopic convention: v2=l2
                    # instead of l2 = [v2::v2+1::2]
                    # It follows than gvib = v2+1
                    # This is added later
                    for v3 in range(v3max + 1):
                        """
                        Spectroscopic rule
                        ------------------
                        J>=l2: as in the symmetric rotor
                        ------------------
                        """
                        viblvl = vib_lvl_name(v1, v2, l2, v3)

                        # First calculate vibrational energy only, for 2-T models
                        if not calc_Evib_harmonic_anharmonic:

                            # recompute
                            # Note Evib = Evib1 + Evib2 + Evib3 should be true
                            # if there is no coupling. General case here: we
                            # recompute each separately
                            if calc_Evib_per_mode:
                                Evib1 = ElecState.Erovib(v1, 0, 0, 0, 0)
                                Evib2 = ElecState.Erovib(0, v2, l2, 0, 0)
                                Evib3 = ElecState.Erovib(0, 0, 0, v3, 0)
                                Evib123 = ElecState.Erovib(v1, v2, l2, v3, 0)
                                # note: the latest could be different from Evib1+Evib2+Evib3
                                # if perturbations are taken into account
                            else:
                                Evib123 = ElecState.Erovib(v1, v2, l2, v3, 0)

                        else:  # calc_Evib_harmonic_anharmonic

                            pass
                            # NotImplemented. Or rather, (G1_h, G1_a), (G2_h, G2_a), (G3_h, G3_a)
                            # is returned anyway by the call to ElecState.Ehaj

                        for J in range(l2, Jmax + 1):
                            #                            roviblvl = '({0}, {1}, {2}, {3}, {4})'.format(v1, v2, l2, v3, J)

                            # Energy level rules

                            # Note: for most isotopes only half the levels exist.
                            # because the other have a gs=0 degeneracy.

                            # Here we still calculate the energy for all
                            # levels because for some reason CDSD features
                            # transitions with some of these levels, and
                            # it creates KeyErrors if we dont
                            # (can there be non null 'forbidden' levels as
                            # there are forbidden transitions?).
                            # Anyway, their linestrengths are very small,
                            # and they wont be acounted for in the
                            # partition function because of gs=0

                            if not calc_Evib_harmonic_anharmonic:

                                # harmonic, corrected
                                Etot = ElecState.Erovib(v1, v2, l2, v3, J)
                                if Etot > Ediss:
                                    break
                                # recompute
                                # Note Evib = Evib1 + Evib2 + Evib3 should be true
                                # if there is no coupling. General case here: we
                                # recompute each separately
                                if calc_Evib_per_mode:
                                    # the following has been precomputed before the
                                    # rotational loop
                                    #                                    Evib1 = ElecState.Erovib(v1, 0, 0, 0, 0)
                                    #                                    Evib2 = ElecState.Erovib(0, v2, l2, 0, 0)
                                    #                                    Evib3 = ElecState.Erovib(0, 0, 0, v3, 0)
                                    #                                    Evib = ElecState.Erovib(v1, v2, l2, v3, 0)

                                    # redefining 'columns' names at each iteration, but
                                    # there is less risk to invert names and data
                                    columns = [
                                        "v1",
                                        "v2",
                                        "l2",
                                        "v3",
                                        "j",
                                        "viblvl",
                                        "E",
                                        "Evib123",
                                        "Evib1",
                                        "Evib2",
                                        "Evib3",
                                    ]
                                    levels.append(
                                        [
                                            v1,
                                            v2,
                                            l2,
                                            v3,
                                            J,
                                            viblvl,
                                            Etot,
                                            Evib123,
                                            Evib1,
                                            Evib2,
                                            Evib3,
                                        ]
                                    )
                                else:
                                    # the following has been precomputed before the
                                    # rotational loop
                                    #                                    Evib = ElecState.Erovib(v1, v2, l2, v3, 0)

                                    # redefining 'columns' names at each iteration, but
                                    # there is less risk to invert names and data
                                    columns = [
                                        "v1",
                                        "v2",
                                        "l2",
                                        "v3",
                                        "j",
                                        "viblvl",
                                        "E",
                                        "Evib123",
                                    ]
                                    levels.append(
                                        [v1, v2, l2, v3, J, viblvl, Etot, Evib123]
                                    )

                            else:  # calc_Evib_harmonic_anharmonic

                                if calc_Evib_per_mode:
                                    (
                                        (G1_h, G1_a),
                                        (G2_h, G2_a),
                                        (G3_h, G3_a),
                                        FJ,
                                    ) = ElecState.Ehaj(v1, v2, l2, v3, J)
                                    Evib1 = G1_h + G1_a
                                    Evib2 = G2_h + G2_a
                                    Evib3 = G3_h + G3_a
                                    Evib123 = Evib1 + Evib2 + Evib3
                                    Etot = Evib123 + FJ
                                    if Etot > Ediss:
                                        break
                                    # recompute
                                    # store

                                    # redefining 'columns' names at each iteration, but
                                    # there is less risk to invert names and data
                                    columns = [
                                        "v1",
                                        "v2",
                                        "l2",
                                        "v3",
                                        "j",
                                        "viblvl",
                                        "E",
                                        "Evib123",
                                        "Evib1",
                                        "Evib2",
                                        "Evib3",
                                        "Evib1_h",
                                        "Evib1_a",
                                        "Evib2_h",
                                        "Evib2_a",
                                        "Evib3_h",
                                        "Evib3_a",
                                    ]
                                    levels.append(
                                        [
                                            v1,
                                            v2,
                                            l2,
                                            v3,
                                            J,
                                            viblvl,
                                            Etot,
                                            Evib123,
                                            Evib1,
                                            Evib2,
                                            Evib3,
                                            G1_h,
                                            G1_a,
                                            G2_h,
                                            G2_a,
                                            G3_h,
                                            G3_a,
                                        ]
                                    )
                                else:
                                    raise NotImplementedError

                            Jmax_calc = max(J, Jmax_calc)
        pb.done()

        df = pd.DataFrame(levels, columns=columns)

        # Calculate missing energies
        # --------------------------

        df["Erot"] = df.E - df.Evib123

        if group_energy_modes_in_2T_model == (["Evib1", "Evib2", "Evib3"], ["Erot"]):
            df["Evib"] = df.Evib123  # + Evib2 + Evib3
            df["Erot"] = df.Erot

        elif group_energy_modes_in_2T_model == (["Evib3"], ["Evib1", "Evib2", "Erot"]):
            if not all_in(["Evib1", "Evib2", "Evib3"], df):
                raise KeyError(
                    "You need Evib1, Evib2, Evib3 calculated separately. "
                    + "Use calc_Evib_per_mode=True"
                )
            df["Evib"] = df.Evib3
            df["Erot"] = df.Evib1 + df.Evib2 + df.Erot

        else:
            raise NotImplementedError(
                "group_energy_modes_in_2T_model: {0}".format(
                    group_energy_modes_in_2T_model
                )
            )
        #         # same as above but harder to read
        #        Evib_columns, Erot_columns = group_energy_modes_in_2T_model
        #        df['Evib20'] = pd.concat([df[k] for k in Evib_columns])
        #        df['Erot2'] = pd.concat([df[k] for k in Erot_columns])
        assert np.allclose(df.Evib + df.Erot, df.E)

        # Get Degeneracies
        # ----------------

        # ... remember: HITRAN convention with v2=l2 -> gvib degeneracy is "v2+1"
        # TODO: make a function call to radis.db.degeneracies.gvib ?
        df["gvib"] = df.v2 + 1
        df["gj"] = 2 * df.j + 1

        gi = self.gi(ElecState)
        gs = self.gs(ElecState)

        # %%
        def is_symmetric(v1, v2, l2, v3):
            """Returns whether a CO2 ``v1v2'l2'v3`` vibrational level is
            symmetric.

            Notes
            -----

            Spectroscopic rules for CO2:

            - Because ground state CO2 is symmetric, and oxygen ^{16}O
             has no spin, for the CO2 1st isotope (626) the negative
             rotational levels* are missing

            *We remind than rotational levels are called positive/negative,
            if the total eigenfunction (including all modes) is
            unchanged/changed, respectively, by reflection of all nuclei
            and electrons at the origin.

            Hence:

            - symmetric vibrational levels -> only even j numbers
            - assymmetric vibrational levels -> only odd j numbers

            If oxygen atoms are different isotopologues (ex: CO2 628), then all
            levels exist

            References
            ----------

            See Section (2.3) in The CO2 Laser by Witteman 1987,
            ISBN 3540477446, 9783540477440

            Exemples
            --------

            Output of is_symmetric(v1, v2, l2, v3) for different levels::

                >>> is_symmetric(1,0,0,0)
                True
                >>> is_symmetric(2,0,0,0)
                True
                >>> is_symmetric(0,2,0,0)
                True
                >>> is_symmetric(0,2,2,0)
                True
                >>> is_symmetric(0,0,0,2)
                True
                >>> is_symmetric(0,0,0,1)
                False
                >>> is_symmetric(0,2,0,1)
                False
                >>> is_symmetric(1,0,0,1)
                False
                >>> is_symmetric(0,3,1,0)
                False
            """
            sym = (-1) ** v2 * (-1) ** v3
            return sym == 1

        def is_even(J):
            """Return whether J is an even or odd rotational level."""
            #            return not bool(J % 2)   # works for float only
            return 1 - np.mod(J, 2)

        if isinstance(gs, tuple):
            # Different degeneracy for symmetric and antisymmetric rotational
            # levels for this isotope. We need the symmetry of rotational levels!

            # Get symmetry of rotational levels

            # Spectroscopic rules
            # -------------------
            # for a spin-symmetric isotope (CO2 626):
            # - if vibrational level is asymmetric, even (+) rotational
            #   levels do not exist
            # - if vibrational level is symmetric, odd (-) rotational
            #   levels do not exist
            # -------------------

            # Get vibrational symmetry: True if Symmetric, False if Antisymmetric
            viblvl_is_symmetric = is_symmetric(
                df.v1.values, df.v2.values, df.l2.values, df.v3.values
            )
            # Get rotational parity: True if even, False if odd
            rotlvl_is_even = is_even(df.j.values)

            # Derive rotational symmetry: True if Symmetric
            rotlvl_is_symmetric = viblvl_is_symmetric == rotlvl_is_even

            # Get state-dependant parity for symmetric and antisymmetric levels
            gs_s, gs_a = gs
            gs = np.ones_like(df.j.values) * gs_a
            gs[rotlvl_is_symmetric] = gs_s

        # Calculate rotational degeneracy
        df["grot"] = gs * gi * df.gj

        #        # just checking:
        #        if not (df.Evib == (df.Evib1 + df.Evib2 + df.Evib3)).all():
        #            raise ValueError('Vib energy is not the sum of all mode energies')
        # not necessarily the case anymore: could be Evib, Erot = Evib3, Evib1+Evib2+Erot

        if self.verbose:
            print(
                "Database generated up to v1={0}, v2={1}, v3={2}, J={3}".format(
                    v1max, v2max, v3max, Jmax_calc
                )
            )

        self.df = df

    def gs(self, ElecState):  # , viblvl):
        """Get state specific rotational degeneracy.

        Parameters
        ----------

        ElecState
            an ElectronicState, that contains molecule id and isotope number

        See Also
        --------

        :func:`~radis.db.degeneracies.gs`
        """

        from radis.db.degeneracies import gs

        M, I = ElecState.id, ElecState.iso
        return gs(M, I)

    def gi(self, ElecState):
        """Get state independant rotational degeneracy. Typically depends on
        the isotope.

        See Also
        --------

        :func:`~radis.db.degeneracies.gi`
        """
        from radis.db.degeneracies import gi

        M, I = ElecState.id, ElecState.iso
        return gi(M, I)


# %% Test
if __name__ == "__main__":

    from radis.test.levels.test_partfunc import _run_testcases

    print("Testing parfunc: {0}".format(_run_testcases()))
