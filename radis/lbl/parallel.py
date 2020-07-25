# -*- coding: utf-8 -*-
"""
Created on Mon May 22 18:23:18 2017

@author: erwan

Build a Parallel Spectrum Factory

Summary
-------

These loads the database once, then calculates a given number of cases in 
parallel using a dedicated SpectrumFactory. Each case is calculated on one 
core only. This is more suited to calculations of large number of small spectra. 
If you want to calculate only one large spectra consider using a single SpectrumFactory
in parallel mode

Routine Listing
---------------

Most methods are written in inherited class with the following inheritance scheme:
    
:py:class:`~radis.lbl.loader.DatabankLoader` > :py:class:`~radis.lbl.base.BaseFactory` > 
:py:class:`~radis.lbl.broadening.BroadenFactory` > :py:class:`~radis.lbl.bands.BandFactory` > 
:py:class:`~radis.lbl.factory.SpectrumFactory` > :py:class:`~radis.lbl.parallel.ParallelFactory`

.. inheritance-diagram:: radis.lbl.parallel.ParallelFactory
   :parts: 1

Examples
--------

Just replace SpectrumFactory by ParallelFactory in your codes. And have a list
as an input of `eq_spectrum` or `non_eq_spectrum`. 

See Also
--------

Refer to :class:`~radis.lbl.factory.SpectrumFactory` for more information. 

-------------------------------------------------------------------------------

"""

from __future__ import absolute_import, print_function, unicode_literals, division
from radis.lbl import SpectrumFactory
from radis.misc.basics import is_list, is_float
from radis.misc.utils import FileNotFoundError
from multiprocessing import Pool, cpu_count

# from multiprocessing.pool import ThreadPool
import sys
from time import time
import numpy as np
from copy import deepcopy
from warnings import warn
import os
from six.moves import zip
from uuid import uuid1
from radis.lbl.loader import df_metadata
from radis.misc.basics import expand_metadata


def _distribute_eq_spectrum(args):
    """ Clone the Factory and calculate eq_spectrum over several clones """
    cast_factory, Tgas, mole_fraction, path_length = args
    # ... (dev) must match p.map(_distribute_eq_spectrum...) order
    factory = deepcopy(cast_factory)
    # Update id
    factory._id = uuid1()
    return SpectrumFactory.eq_spectrum(
        factory, Tgas, mole_fraction=mole_fraction, path_length=path_length
    )


def _distribute_noneq_spectrum(args):
    """ Clone the Factory and calculate non_eq_spectrum over several clones """
    cast_factory, Tvib, Trot, mole_fraction, path_length = args
    # ... (dev) must match p.map(_distribute_noneq_spectrum...) order
    factory = deepcopy(cast_factory)
    # Update id
    factory._id = uuid1()
    return SpectrumFactory.non_eq_spectrum(
        factory, Tvib, Trot, mole_fraction=mole_fraction, path_length=path_length
    )


class ParallelFactory(SpectrumFactory):
    """ A Parallel version of :class:`~radis.lbl.parallel.SpectrumFactory`. 
    Use `Nprocs` as argument to change the number of processors to be used.


    Notes
    -----

    for Windows users:

    On Windows the subprocesses will import (i.e. execute) the main module at start. 
    You need to protect the main code like this to avoid creating subprocesses recursively::

        from radis.lbl import ParallelFactory
        if __name__ == '__main__':    
            # your code

    More info: https://stackoverflow.com/questions/18204782/runtimeerror-on-windows-trying-python-multiprocessing

    Parallel is only implemented on Python 3. 

    Parameters
    ----------

    wavenum_min: (cm-1)
        minimum wavenumber to be processed in cm^-1

    wavenum_max: (cm-1)
        maximum wavenumber to be processed in cm^-1

    Tref: (K)
        reference temperature for calculations, HITRAN database uses 296 Kelvin
        default: 3400K

    pressure: (bar)
        partial pressure of gas in bar. Default 1.01325 (1 atm)

    mole_fraction: N/D
        species mole fraction. Default 1. Note that the rest of the gas
        is considered to be air for collisional broadening.

    path_length: (cm)
        path length in cm. Default 1.

    isotope: int, str, or list of int
        isotope id (sorted by relative density). Default [1,2] (eg: CO2-626, CO2-636 for CO2)

    Other Parameters
    ----------------

    Nprocs: int
        Number of processors to use. Default total number of threads. 

    Tref: K
        Reference temperature for calculations (linestrength temperature
        correction). HITRAN database uses 296 Kelvin. Default 296 K

    broadening_max_width: cm-1
        Full width over which to compute the broadening. Large values will create
        a huge performance drop (because convolutions are not vectorized).
        Also, the calculated spectral range is increased (by broadening_max_width/2
        on each side) to take into account overlaps from out-of-range lines. 
        Default 10 cm-1.

    wstep: cm-1
        Spacing of calculated absorbance spectrum. Default 0.01 cm-1

    cutoff: float (~ unit of Linestrength: cm-1/(#.cm-2))
        discard linestrengths that are lower that this, to reduce calculation
        times. 1e-27 is what is generally used to generate databases such as 
        CDSD. If 0, no cutoff. Default 1e-27.

    bplot: boolean
        plot intermediary results (like slit function generation). Default ``False``.

    save_memory: boolean
        if ``True``, removes databases calculated by intermediate functions (for
        instance, delete the full database once the linestrength cutoff criteria
        was applied). This saves some memory but requires to reload the database
        & recalculate the linestrength for each new parameter. Default ``False``.

    chunksize: float
        Splits the lines database in several chuncks during calculation, else
        the multiplication of lines over all spectral range takes too much memory
        and slows the system down. Chuncksize let you change the default chunck
        size. Default 3e7
        
    Examples
    --------
    
    Refer to the online :ref:`Examples <label_examples>`. 


    See Also
    --------

    :class:`~radis.lbl.factory.SpectrumFactory`
    """

    def __init__(self, *args, **kwargs):

        Nprocs = kwargs.pop("Nprocs", cpu_count())  # default to cpu_count()
        kwargs["parallel"] = False  # calculate each spectrum on one core only

        super(ParallelFactory, self).__init__(*args, **kwargs)

        self.misc.Nprocs = Nprocs

        if os.name == "nt" and self.verbose:
            print(
                (
                    "Warning. On Windows, if executing the parallel versions of "
                    + "eq_spectrum and non_eq_spectrum in a script, make sure you "
                    + "embed them in `if __name__ == '__main__'` or daemons processes"
                    + " won't start"
                )
            )

    def eq_spectrum(self, Tgas, mole_fraction=None, path_length=None):
        """ Generate a spectrum at equilibrium

        Parameters
        ----------

        Tgas: list or float
            Gas temperature (K)

        mole_fraction: list or float
            database species mole fraction. If None, Factory mole fraction is used.

        path_length: list or float
            slab size (cm). If None, Factory mole fraction is used.

        Returns
        -------

        Returns a :class:`~radis.spectrum.spectrum.Spectrum` object

        Use the :meth:`~radis.spectrum.spectrum.Spectrum.get` method to get something
        among ``['radiance', 'radiance_noslit', 'absorbance', etc...]``

        Or directly the :meth:`~radis.spectrum.spectrum.Spectrum.plot` method
        to plot it

        See [1]_ to get an overview of all Spectrum methods

        References
        ----------

        .. [1] RADIS doc: `Spectrum how to? <https://radis.readthedocs.io/en/latest/spectrum/spectrum.html#label-spectrum>`__


        Notes
        -----

        Process:

        Calculate line strenghts correcting the CDSD reference one. Then call
        the main routine that sums over all lines


        """

        args = self._format_input(
            **{"Tgas": Tgas, "mole_fraction": mole_fraction, "path_length": path_length}
        )
        Tgas = args.pop("Tgas")
        mole_fraction = args.pop("mole_fraction")
        path_length = args.pop("path_length")

        # Expand metadata before copy, else it will be lost.
        # @dev: See :py:data:`~radis.lbl.loader.df_metadata` for more explanations
        expand_metadata(self.df0, [k for k in df_metadata if hasattr(self.df0, k)])

        N = len(Tgas)
        cast_factory = [self] * N  # note that this is the copy of a same object!

        Nprocs = min(self.misc.Nprocs, N)
        if self.verbose:
            print(("Calculate {0} spectra on {1} procs".format(len(Tgas), Nprocs)))
        t1 = time()

        try:
            with Pool(Nprocs) as p:  # Python 3 only
                spec_list = p.map(
                    _distribute_eq_spectrum,
                    list(zip(cast_factory, Tgas, mole_fraction, path_length)),
                )
        except AttributeError:
            if sys.version_info[0] == 2:
                raise NotImplementedError(
                    "parallel is only implemented in Python 3. Time to switch!"
                )
            else:
                raise

        if self.verbose:
            print(("{0} spectra calculated in {1:.1f}s".format(N, time() - t1)))

        return spec_list

    def non_eq_spectrum(self, Tvib, Trot, mole_fraction=None, path_length=None):
        """ Calculate emission spectrum in non-equilibrium case. Calculates
        absorption with broadened linestrength and emission with broadened
        Einstein coefficient.

        Parameters
        ----------

        Tvib: list or float
            vibrational temperature [K]

        Trot: list or float
            rotational temperature [K]

        mole_fraction: list or float
            database species mole fraction. If None, Factory mole fraction is used.

        path_length: list or float. If None, Factory mole fraction is used.
            slab size (cm)

        Returns
        -------

        Returns a :class:`~radis.spectrum.spectrum.Spectrum` object

        Use the :meth:`~radis.spectrum.spectrum.Spectrum.get` method to get something
        among ``['radiance', 'radiance_noslit', 'absorbance', etc...]``

        Or directly the :meth:`~radis.spectrum.spectrum.Spectrum.plot` method
        to plot it

        See [1]_ to get an overview of all Spectrum methods

        References
        ----------

        .. [1] RADIS doc: `Spectrum how to? <https://radis.readthedocs.io/en/latest/spectrum/spectrum.html#label-spectrum>`__


        Notes
        -----

        Hypothesis:

        - Trot = Ttrans

        """
        args = self._format_input(
            **{
                "Tvib": Tvib,
                "Trot": Trot,
                "mole_fraction": mole_fraction,
                "path_length": path_length,
            }
        )
        Tvib = args["Tvib"]
        Trot = args["Trot"]
        mole_fraction = args["mole_fraction"]
        path_length = args["path_length"]

        # Expand metadata before copy, else it will be lost.
        # @dev: See :py:data:`~radis.lbl.loader.df_metadata` for more explanations
        expand_metadata(self.df0, [k for k in df_metadata if hasattr(self.df0, k)])

        N = len(Tvib)
        cast_factory = [self] * N  # note that this is the copy of a same object!

        Nprocs = min(self.misc.Nprocs, N)
        if self.verbose:
            print(("Calculate {0} spectra on {1} procs".format(len(Tvib), Nprocs)))
        t1 = time()

        try:
            with Pool(Nprocs) as p:  # Note: Python 3 syntax only
                spec_list = p.map(
                    _distribute_noneq_spectrum,
                    list(zip(cast_factory, Tvib, Trot, mole_fraction, path_length)),
                )
        except AttributeError:
            if sys.version_info[0] == 2:
                raise NotImplementedError(
                    "parallel is only implemented in Python 3. Time to switch!"
                )
            else:
                raise

        if self.verbose:
            print(("{0} spectra calculated in {1:.1f}s".format(N, time() - t1)))

        return spec_list

    def _format_input(self, **kwargs):
        """ Format the input of parallel calculation request. Returns lists of 
        same lengths that can be parsed with zip(). User input can be lists, 
        or floats instead of constant-values list """
        N = None  # list length if there are list involved
        kwout = {}
        for k, v in kwargs.items():
            if is_list(v):
                if N is None:
                    N = len(v)
                else:
                    if len(v) != N:
                        raise ValueError(
                            "Invalid input for {0}: ".format(k)
                            + "all list should have the same length"
                            "(you may use floats too)"
                        )
                kwout[k] = v

        # Now let's turn floats into list of length N
        if N is None:
            warn("Using ParallelFactory for single case")
            N = 1
        for k, v in kwargs.items():
            if type(v) in [int, np.int64]:
                v = float(v)  # int is not serializable (for some reason)
            if is_list(v):
                continue  # already done
            elif is_float(v) or v is None:
                kwout[k] = [v] * N
            else:
                raise ValueError(
                    "Invalid input for {0}: ".format(k)
                    + "input should be list-like or float-like"
                    + "({0} instead)".format(type(v))
                )

        # check output is correct:
        for v in kwout.values():
            assert len(v) == N

        return kwout


# %% Main
if __name__ == "__main__":

    from radis.test.lbl.test_parallel import _run_testcases

    print(("Tested parallel.py:", _run_testcases(plot=True)))
