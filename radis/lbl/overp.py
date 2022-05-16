# -*- coding: utf-8 -*-
"""Created on Thu Sep 14 13:44:35 2017.

@author: erwan

Summary
-------

An experimental module to merge precalculated vibrational bands in post-processing,
and recompute a new spectrum at a new temperature or with overpopulations, without
having to recalculate the broadening of each line

.. warning::

    Order of magnitude faster, but only valid under optically thin conditions as the
    rescaling of absorption doesnt scale induced emission properly

----------
"""

from time import time
from warnings import warn

import numpy as np
import pandas as pd
from numpy import exp

try:  # Proper import
    from .labels import vib_lvl_name_cdsd_pc, vib_lvl_name_cdsd_pcN
except ImportError:  # if ran from here
    from radis.lbl.labels import vib_lvl_name_cdsd_pc, vib_lvl_name_cdsd_pcN
from radis.misc.basics import is_float
from radis.phys.constants import hc_k
from radis.spectrum.rescale import (
    rescale_abscoeff,
    rescale_absorbance,
    rescale_emisscoeff,
    rescale_radiance_noslit,
    rescale_transmittance_noslit,
)
from radis.spectrum.spectrum import Spectrum

# keys that may not be equals for different bands
_IGNORE_KEYS = ["band", "band_htrn", "viblvl_u", "viblvl_l"]


class LevelsList(object):
    """A class to generate a Spectrum from a list of precalculated bands at a
    given reference temperature.

    .. warning::

        only valid under optically thin conditions!!

    See Also
    --------
    """

    # hardcode attribute names, but can save a lot of memory
    __slots__ = [
        "Tvib_ref",
        "Trot_ref",
        "vib_distribution_ref",
        "mole_fraction_ref",
        "path_length_ref",
        "sorted_bands",
        "cum_weight",
        "E_bands",
        "bands_ref",
        "verbose",
        "bands",
        "parfunc",
        "levelsfmt",
        "lvl_index",
        "vib_levels",
        "copy_lines",
        "Qref",
        "Qvib_ref",
    ]

    def __init__(
        self, parfunc, bands, levelsfmt, sortby="Ei", copy_lines=False, verbose=True
    ):
        """
        Parameters
        ----------
        bands: dict of bands
            bands are Spectrum objects calculated at equilibrium or non-equilibrim.
            Tgas, or (Tvib, Trot) must be given and the same in all bands conditions.

        """

        # Check inputs
        for br, s in bands.items():
            if not isinstance(s, Spectrum):
                raise ValueError(
                    "`bands` must be a list of Spectrum objects. "
                    + "Got {0}".format(type(s))
                )
            if s.lines is None:
                raise ValueError(
                    "To be used in LevelsList spectra must have been calculated with 'export_lines=True'"
                )

        # Assert all conditions are the same
        bands_list = list(bands.values())
        cond_ref = bands_list[0].conditions
        for s in bands_list[1:]:  # type(s): Spectrum
            if list(s.conditions.keys()) != list(cond_ref.keys()):
                raise ValueError("Input conditions must be the same for all bands")
            for k in cond_ref.keys():
                if k in _IGNORE_KEYS:  # don't compare these entries
                    continue
                if s.conditions[k] != cond_ref[k]:
                    raise ValueError(
                        "Input conditions must be the same for all bands"
                        + ". Got a difference for entry {0}: {1} vs {2}".format(
                            k, s.conditions[k], cond_ref[k]
                        )
                    )

        if "Tvib" in cond_ref:
            Tvib_ref = cond_ref["Tvib"]
        else:
            Tvib_ref = cond_ref["Tgas"]
        if "Trot" in cond_ref:
            Trot_ref = cond_ref["Trot"]
        else:
            Trot_ref = cond_ref["Tgas"]
        vib_distribution_ref = cond_ref["vib_distribution"]
        if "overpopulation" in cond_ref:
            raise NotImplementedError

        # stores input conditions
        self.Tvib_ref = Tvib_ref
        self.Trot_ref = Trot_ref
        self.mole_fraction_ref = cond_ref["mole_fraction"]
        self.path_length_ref = cond_ref["path_length"]
        self.vib_distribution_ref = vib_distribution_ref

        # stores computation parameters
        # TODO: store it in parfunc instead (and make parfunc an EnergyDatabase)
        self.levelsfmt = levelsfmt

        # %% Misc
        self.verbose = verbose
        self.copy_lines = copy_lines

        # %% Set up energies  ( use a partition function object for the moment)

        self.parfunc = parfunc
        self._init_levels(sortby)
        self._connect_bands_to_levels(bands)

        Qvib = self._calc_vib_populations(
            Tvib=self.Tvib_ref, vib_distribution=self.vib_distribution_ref
        )

        self.Qref = parfunc.at_noneq(
            Tvib_ref, Trot_ref, vib_distribution=vib_distribution_ref
        )  # TODO: add overpopulations in partition function
        self.Qvib_ref = Qvib

        self.vib_levels["nvib_ref"] = self.vib_levels["nvib"]
        del self.vib_levels["nvib"]

        # %% Set up bands
        self._init_bands_ref(bands, sortby)

    def _init_bands_ref(self, bands, sortby):

        # %% Parse bands

        #        # Get band names   (if bands was a list)
        #        band_dict = {}
        #        for bd in bands:
        #            try:
        #                band_name = set(bd.lines['band'])
        #            except:
        #                raise ValueError('Bands must have been created with a `band` attribute'+\
        #                                 ' in lines')
        #            assert len(band_name) == 0
        #            band_dict[band_name[0]] = bd
        #        self.bands = band_dict
        #        self.bands0 = {br: bands[br].copy(copy_lines=False) for br in bands}  # make a copy

        # Get energy of upper level of band
        E_bands = {}
        for br, band in bands.items():  # type(band): Spectrum
            if br == "others":
                continue
            E_bands[br] = np.mean(list(set(band.lines["Evibu"])))
            # Todo: should be all equals ?
            assert len(set(band.lines["Evibu"])) == 1

        # Calculate band weight
        weight = {}
        for br, band in bands.items():  # type(band): Spectrum
            weight[br] = band.lines[sortby].sum()
        tot = np.sum([np.array(w) for w in weight.values()])
        for br, band in bands.items():
            weight[br] /= tot
        sorted_bands = sorted(bands, key=lambda br: weight[br])[::-1]
        cum_weight = np.array([weight[br] for br in sorted_bands])
        cum_weight = np.cumsum(cum_weight)

        # Store
        self.sorted_bands = sorted_bands
        self.cum_weight = cum_weight
        self.E_bands = E_bands

        #        # Delete lines (makes copies faster)
        #        for br in self.bands0:
        #            del self.bands0[br].lines

        # Initialize bands
        self.bands_ref = {
            br: band.copy(copy_lines=self.copy_lines) for br, band in bands.items()
        }

    def _init_levels(self, sortby):
        """
        Notes
        -----

        vib level in (p,c,j,N) notation is ambiguous. Here we use (p,c) based
        on energy difference and perturbation rules
        """

        df = self.parfunc.df
        levelsfmt = self.levelsfmt

        if levelsfmt == "cdsd-pc":
            vib_index = ["p", "c"]

            def vib_lvl_name_cdsd(p, c, N):
                return vib_lvl_name_cdsd_pc(p, c)

        elif levelsfmt == "cdsd-pcN":
            vib_index = ["p", "c", "N"]
            vib_lvl_name_cdsd = vib_lvl_name_cdsd_pcN

        else:
            raise NotImplementedError(levelsfmt)

        vib_levels = df.drop_duplicates(subset=vib_index)  # 192 ms ± 1.13 ms

        # Add levels
        def add_viblvl(row):
            """ """
            #    row['vib_lvl'] = '({p},{c},{N})'.format(**dict([(k,int(row[k])) for k in ['p', 'c', 'N']]))
            #            row['viblvl'] = '({p},{c})'.format(**dict([(k,int(row[k])) for k in ['p', 'c']]))
            row["viblvl"] = vib_lvl_name_cdsd(row.p, row.c, row.N)
            return row

        t0 = time()
        vib_levels = vib_levels.apply(add_viblvl, axis=1)
        print("Added viblvl in {0:.2f}s".format(time() - t0))

        # vib_levels = vib_levels.set_index(['p', 'c', 'N'])
        vib_levels = vib_levels.set_index(["viblvl"])

        # %% Generate index of  levels

        lvl_index = dict().fromkeys(list(vib_levels.index))
        for (
            k
        ) in (
            lvl_index
        ):  # Note: do not fill directly in frmokeys or you'll get the same shared list
            lvl_index[k] = {"bands_where_low": [], "bands_where_up": []}

        self.lvl_index = lvl_index
        self.vib_levels = vib_levels

    def _connect_bands_to_levels(self, s_bands):

        index = self.lvl_index

        for s in s_bands.values():
            if s.name == "others":
                continue

            viblvl_u = s.lines["viblvl_u"].iloc[0]
            viblvl_l = s.lines["viblvl_l"].iloc[0]

            index[viblvl_u]["bands_where_up"].append(s)
            index[viblvl_l]["bands_where_low"].append(s)

    def plot_vib_populations(self, nfig=None, **kwargs):
        """Plot current distribution of vibrational levels.

        By constructions populations are shown as divided by the state degeneracy,
        i.e,            g =   gv   *   (2J+1)   *    gi    *     gs

        Parameters
        ----------

        nfig: str, int
            name of Figure to plot on

        kwargs: **dict
            arguments are forwarded to plot()
        """

        import matplotlib.pyplot as plt

        vib_levels = self.vib_levels

        plt.figure(num=nfig)
        plt.plot(vib_levels.E, vib_levels.nvib, "ok", **kwargs)
        plt.xlabel("Energy (cm-1)")
        plt.ylabel("Population")
        plt.yscale("log")
        plt.xlim(xmin=0)

        return

    def _calc_vib_populations(
        self, Tvib, vib_distribution="boltzmann", overpopulation=None
    ):
        """Calculate vibrational populations from Tref to new Tvib and store
        results in vib_levels dataframe This does not modify the spectra yet!

        By constructions populations are calculated divided by the state degeneracy,
        i.e,            g =   gv   *   (2J+1)   *    gi    *     gs
        This means we should only use ratios for rescaling

        The information on state degeneracy and isotopic abundance is already
        included in the linestrength / emission integral, hence in the pre-calculated
        emisscoeff / abscoeff
        """

        vib_levels = self.vib_levels
        if overpopulation is None:
            overpopulation = {}

        if is_float(Tvib):

            # Get new population
            E_vib = vib_levels["Evib"]
            g = 1  # explicitely calculate populations divided by degeneracy
            # this means we should only use ratios for rescaling
            if vib_distribution == "boltzmann":
                nvibQvib = g * exp(-hc_k * E_vib / Tvib)
            else:
                raise NotImplementedError(
                    "vib_distribution: {0}".format(vib_distribution)
                )

        else:
            Tvib1, Tvib2, Tvib3 = Tvib

            # Get new population
            E_vib1 = vib_levels["Evib1"]
            E_vib2 = vib_levels["Evib2"]
            E_vib3 = vib_levels["Evib3"]
            g = 1  # explicitely calculate populations divided by degeneracy
            # this means we should only use ratios for rescaling
            if vib_distribution == "boltzmann":
                nvibQvib = (
                    g
                    * exp(-hc_k * E_vib1 / Tvib1)
                    * exp(-hc_k * E_vib2 / Tvib2)
                    * exp(-hc_k * E_vib3 / Tvib3)
                )
            else:
                raise NotImplementedError(
                    "vib_distribution: {0}".format(vib_distribution)
                )

        # Add overpopulation
        if overpopulation != {}:
            for k, ov in overpopulation.items():
                nvibQvib.loc[k] *= ov
                # TODO: test
            warn(
                "NotImplemented: partition function overpopulation correction not tested"
            )

        # Normalize with partition function
        Qvib = nvibQvib.sum()
        nvib = nvibQvib / Qvib

        # update dataframe
        vib_levels["nvib"] = nvib
        vib_levels["Qvib"] = Qvib

        return Qvib

    #    def set_cutoff(self, cutoff=1e-3):
    #        ''' Discard most of the bands and keep error below cutoff
    #
    #        Warning: this cannot be undone!'''
    #
    #        cum_weight = self.cum_weight
    #        sorted_bands = self.sorted_bands
    #
    #        cutoff_i = np.argmax((1-cum_weight)<cutoff)
    #        first_bands = sorted_bands[:cutoff_i+1]
    #
    #        # Update bands
    #        new_bands = {br: self.bands_ref[br] for br in first_bands}
    #
    #        if self.verbose:
    #            print(('Set cutoff. Discarded {0}/{1} bands. Estimated error: {2:.2f}%'.format(
    #                len(self.bands_ref)-len(new_bands), len(self.bands_ref),
    #                (1-cum_weight[cutoff_i])*100)))
    #
    #        self.bands_ref = new_bands

    def eq_spectrum(
        self,
        Tgas,
        overpopulation=None,
        mole_fraction=None,
        path_length=None,
        save_rescaled_bands=False,
    ):
        """See :py:meth:`~radis.lbl.factory.SpectrumFactory.eq_spectrum`

        .. warning::

            only valid under optically thin conditions!!

        Parameters
        ----------

        ... same as usually. If None, then the reference value (used to
        calculate bands) is used

        Other Parameters
        ----------------

        save_rescaled_bands: boolean
            save updated bands. Take some time as it requires rescaling all
            bands individually (which is only done on the MergedSlabs usually)
            Default ``False``
        """
        s = self.non_eq_spectrum(
            Tgas,
            Tgas,
            overpopulation=overpopulation,
            mole_fraction=mole_fraction,
            path_length=path_length,
            save_rescaled_bands=save_rescaled_bands,
        )
        del s.conditions["Tvib"]
        del s.conditions["Trot"]
        s.conditions["Tgas"] = Tgas

        return s

    def non_eq_spectrum(
        self,
        Tvib=None,
        Trot=None,
        Ttrans=None,
        vib_distribution="boltzmann",
        overpopulation=None,
        mole_fraction=None,
        path_length=None,
        save_rescaled_bands=False,
    ):
        """See :py:meth:`~radis.lbl.factory.SpectrumFactory.non_eq_spectrum`

        .. warning::

            only valid under optically thin conditions!!

        Parameters
        ----------

        ... same as usually. If None, then the reference value (used to
        calculate bands) is used

        Other Parameters
        ----------------

        save_rescaled_bands: boolean
            save updated bands. Take some time as it requires rescaling all
            bands individually (which is only done on the MergedSlabs usually)
            Default ``False``

        Notes
        -----

        Implementation:

        Generation of a new spectrum is done by recombination of the precalculated
        bands with
        """

        from radis.los import MergeSlabs

        # Restart from a copy each time  (else reference bands are modified by rescaling
        # and we may loose information if rescaling to 0 for instance)
        bands_ref = self.bands_ref
        bands = {br: bands_ref[br].copy(copy_lines=self.copy_lines) for br in bands_ref}

        # Initialize inputs
        if Trot is None:
            Trot = self.Trot_ref
        if Tvib is None:
            Tvib = self.Tvib_ref
        if Trot != self.Trot_ref:
            raise ValueError(
                "Trot {0} doesnt match the reference Trot: {1}".format(
                    Trot, self.Trot_ref
                )
            )
        if mole_fraction is None:
            mole_fraction = self.mole_fraction_ref
        if path_length is None:
            path_length = self.path_length_ref
        if overpopulation is None:
            overpopulation = {}

        #        # Fill missing overpopulation with 1
        #        if overpopulation is None:
        #            overpopulation = dict.fromkeys(list(bands.keys()), 1)
        #        else:
        #            for br in bands.keys():
        #                if not br in overpopulation:
        #                    overpopulation[br] = 1

        #        Tvib_ref = self.Tvib_ref
        #        pop_correction = dict.fromkeys(list(bands.keys()), 1)
        #        if Tvib != Tvib_ref:
        #            raise NotImplementedError
        #            E_bands = self.E_bands
        #            # Correct for vibrational temperature
        #            for br in bands:
        #                pop_correction[br] = exp(-E_bands[br]*hc_k/Tvib)/exp(-E_bands[br]*hc_k/Tvib_ref)

        # Recalculate populations from reference everytime
        vib_levels = self.vib_levels

        # Recalculate partition function
        Qvib = self._calc_vib_populations(
            Tvib=Tvib, vib_distribution=vib_distribution, overpopulation=overpopulation
        )
        # TODO: note: Qvib doesnt take gvib into account. Wrong Qvib but correct populations
        # if rescaled with ratio???
        # ... warning: update_populations set to False here not updated (because we dont
        # ... care much about all levels here. only the one from vib_levels matter)
        if overpopulation is not None:
            Q, _, _ = self.parfunc.at_noneq(
                Tvib,
                Trot,
                overpopulation=overpopulation,
                returnQvibQrot=True,
                update_populations=False,
            )
        else:
            Q = self.parfunc.at_noneq(Tvib, Trot, update_populations=False)

        Qref = self.Qref
        Qvib_ref = self.Qvib_ref

        for br, band in bands.items():  # type(band): Spectrum
            if br == "others":
                continue
            viblvl_u = band.conditions["viblvl_u"]
            viblvl_l = band.conditions["viblvl_l"]
            nu_vib = vib_levels.loc[viblvl_u, "nvib"]
            nu_vib_old = vib_levels.loc[viblvl_u, "nvib_ref"]
            nl_vib = vib_levels.loc[viblvl_l, "nvib"]
            nl_vib_old = vib_levels.loc[viblvl_l, "nvib_ref"]

            # with new, old = (2, 1):
            # n2/n1 = exp(-Evib2/kT2)/exp(-Evib1/kT1) * Qvib1/Qvib2
            #       = n2vib / n1vib * Q1/Q2 * Q2vib / Q1vib
            if vib_distribution == "boltzmann":
                corfactor_u = nu_vib / nu_vib_old * Qref / Q * Qvib / Qvib_ref
                corfactor_l = nl_vib / nl_vib_old * Qref / Q * Qvib / Qvib_ref
            else:
                raise NotImplementedError(
                    "vib_distribution: {0}".format(vib_distribution)
                )

            rescale_updown_levels(band, corfactor_u, 1, corfactor_l, 1)

            # Update lines
            if self.copy_lines:
                band.lines["nu"] *= corfactor_u
                band.lines["nl"] *= corfactor_l
                band.lines["Qvib"] = Qvib
                band.lines["nu_vib"] *= corfactor_u
                band.lines["nl_vib"] *= corfactor_l
                band.lines["Ei"] *= np.nan
                band.lines["S"] = np.nan

        # Get total spectrum
        s = MergeSlabs(*list(bands.values()))

        # populations
        s.populations = pd.DataFrame(vib_levels[["nvib", "Evib"]])

        # rebuild lines
        if self.copy_lines:
            s.lines = pd.concat([band.lines for band in bands.values()])

        # Update total mole fraction if it changed
        mole_fraction_ref = self.mole_fraction_ref
        if mole_fraction != mole_fraction_ref:
            s.rescale_mole_fraction(
                mole_fraction,
                mole_fraction_ref,
                ignore_warnings=True,  # mole_fraction is 'N/A'
                force=True,
            )

        path_length_ref = self.path_length_ref
        if path_length != path_length_ref:
            s.rescale_path_length(path_length, path_length_ref, force=True)

        # Add parameters in conditions:
        s.conditions["overpopulation"] = overpopulation
        s.conditions["mole_fraction"] = mole_fraction  # was 'N/A' after Merge
        # because it's different for all bands
        s.conditions["path_length"] = path_length
        s.conditions["Tvib"] = Tvib
        s.conditions["Trot"] = Trot

        if save_rescaled_bands:
            for br, band in bands.items():  # type(band): Spectrum
                if mole_fraction != mole_fraction_ref:
                    band.rescale_mole_fraction(
                        mole_fraction,
                        mole_fraction_ref,
                        ignore_warnings=True,  # mole_fraction is 'N/A'
                        force=True,
                    )
                if path_length != path_length_ref:
                    band.rescale_path_length(path_length, path_length_ref, force=True)
            self.bands = bands

        return s


# %% Overpopulation functions


def rescale_updown_levels(
    spec, new_nu, old_nu, new_nl, old_nl, ignore_warnings=False, force=False
):
    """Update spectrum with new molar fraction for upper and lower levels.

    Convoluted values (with slit) are dropped in the process.

    Rescales with a ratio new_nu/old_nu.

    This is only valid for emission quantities, under optically thin conditions,
    as rescaling doesnt correct for induced emission.

    .. warning::
        experimental feature

    Parameters
    ----------
    new_nu: float
        new upper state mole fraction
    old_nu: float
        current upper state mole fraction
    new_nl: float
        new lower state mole fraction
    old_nl: float
        current lower state mole fraction

    Other Parameters
    ----------------
    force: boolean
        if False, won't allow rescaling to 0 (not to loose information).
        Default ``False``

    Notes
    -----

    Implementation:

    similar to rescale_path_length() but we have to scale abscoeff & emisscoeff
    Note that this is valid only for small changes in mole fractions. Then,
    the change in line broadening becomes significant and the whole band
    should be recomputed

    IMPORTANT: if editing make sure you use the proper nu and nl. In particular
    when infering emission quantities from absorption quantities this may
    ends up in error in overpopulation rescaling.
    """
    # TODO: Add warning when too large rescaling

    # Check inputs
    # ---------
    #        if old_mole_fraction is not None:
    #            try:
    #                if spec.conditions['mole_fraction'] != old_mole_fraction and not ignore_warnings:
    #                    warn('mole_fraction ({0}) doesnt match value given in conditions ({1})'.format(
    #                            old_mole_fraction, spec.conditions['mole_fraction']))
    #            except KeyError: # mole fraction not defined
    #                pass

    #        else:
    #            try:
    #                old_mole_fraction = spec.conditions['mole_fraction']
    #            except KeyError:
    #                raise KeyError('mole_fraction has to be defined in conditions (or use'+\
    #                                ' `from_mole_fraction`)'),

    if (new_nu < 0 or new_nl < 0) and not force:
        raise ValueError("mole_fraction cannot be negative")
    if (new_nu == 0 or new_nl == 0) and not force:
        raise ValueError(
            "Rescaling to 0 will loose information. Choose force " "= True"
        )

    for q in ["transmittance", "radiance"]:
        qns = q + "_noslit"
        qties = spec.get_vars()
        if q in qties and qns not in qties and not force:
            raise KeyError(
                "Cant rescale {0} if {1} not stored".format(q, qns)
                + " Use force=True to rescale anyway. {0}".format(q)
                + " will be deleted"
            )

    # Get path length
    if "path_length" in list(spec.conditions.keys()):
        path_length = spec.conditions["path_length"]
        true_path_length = True
    else:
        path_length = 1
        true_path_length = False

    # %%  Rescale

    optically_thin = spec.is_optically_thin()

    # Choose which values to recompute
    # ----------
    initial = spec.get_vars()  # quantities initialy in spectrum
    recompute = list(initial)  # quantities to recompute
    rescaled = {}  # quantities rescaled

    if "radiance_noslit" in initial and not optically_thin:
        recompute.append("emisscoeff")
    if (
        "absorbance" in initial
        or "transmittance_noslit" in initial
        or "radiance_noslit" in initial
        and not optically_thin
    ):
        recompute.append("abscoeff")
    extra = []
    recompute = set(recompute)  # remove duplicates

    # Get units
    units = spec.units.copy()

    # Recompute!
    # ----------
    waveunit = spec.get_waveunit()  # keep all quantities in same waveunit

    if "abscoeff" in recompute:
        rescaled, units = rescale_abscoeff(
            spec,
            rescaled,
            initial,  # TODO: remove rescaled = ... Dict is mutable
            old_nl,
            new_nl,
            path_length,
            waveunit,
            units,
            extra,
            true_path_length,
            assume_equilibrium=False,
        )

    if "emisscoeff" in recompute:
        rescaled, units = rescale_emisscoeff(
            spec,
            rescaled,
            initial,
            old_nu,
            new_nu,
            path_length,
            optically_thin,
            waveunit,
            units,
            extra,
            true_path_length,
        )

    if "absorbance" in recompute:
        rescaled, units = rescale_absorbance(
            spec,
            rescaled,
            initial,
            old_nl,
            new_nl,
            path_length,
            path_length,
            waveunit,
            units,
            extra,
            true_path_length,
        )

    if "transmittance_noslit" in recompute:
        rescaled, units = rescale_transmittance_noslit(
            spec,
            rescaled,
            initial,
            old_nl,
            new_nl,
            path_length,
            path_length,
            waveunit,
            units,
            extra,
            true_path_length,
        )

    if "radiance_noslit" in recompute:
        rescaled, units = rescale_radiance_noslit(
            spec,
            rescaled,
            initial,
            old_nu,
            new_nu,
            path_length,
            path_length,
            optically_thin,
            waveunit,
            units,
            extra,
            true_path_length,
        )

    # Save (only) the ones that were in the spectrum initially
    for q in rescaled:
        if q in initial:
            spec._q[q] = rescaled[q]

    # Update units
    for k, u in units.items():
        spec.units[k] = u

    # Drop convoluted values
    for q in ["transmittance", "radiance"]:
        if q in list(spec._q.keys()):
            del spec._q[q]

    # Reapply slit if possible
    if (
        "slit_function" in spec.conditions
        and "slit_unit" in spec.conditions
        and "norm_by" in spec.conditions
    ):
        slit_function = spec.conditions["slit_function"]
        slit_unit = spec.conditions["slit_unit"]
        norm_by = spec.conditions["norm_by"]
        try:
            shape = spec.conditions["shape"]
        except KeyError:
            shape = "triangular"
        spec.apply_slit(
            slit_function=slit_function, unit=slit_unit, shape=shape, norm_by=norm_by
        )

    # Update conditions
    spec.conditions["mole_fraction"] = new_nl if new_nl == new_nu else "N/A"


# def _test(verbose=True, plot=False, *args, **kwargs):
#    """ Generate scalable LevelList at 1500K. Rescale at 2000K and compare with
#    spectrum directly calculated at 2000K
#    """
#
#    from radis.lbl import SpectrumFactory
#
#    iso = 1
#    sf = SpectrumFactory(
#        wavelength_min=4170,
#        wavelength_max=4200,
#        mole_fraction=1,
#        path_length=0.05,
#        cutoff=1e-25,
#        #                     isotope=[1,2],
#        isotope=iso,
#        db_use_cached=True,
#        wstep=0.01,
#        broadening_max_width=10,
#        medium="air",
#        verbose=verbose,
#    )
#
#    sf.load_databank("CDSD")
#
#    parfunc = sf.parsum_calc["CO2"][iso]["X"]
#
#    # %% Fill levels for all bands
#
#    Tref = 1500
#    s_bands = sf.non_eq_bands(Tvib=Tref, Trot=Tref)
#    lvlist = LevelsList(parfunc, s_bands, sf.params.levelsfmt)
#
#    # %% Test
#    Tvib = 2000
#    s_resc = lvlist.non_eq_spectrum(Tvib=Tvib, Trot=Tref)
#    s0 = sf.non_eq_spectrum(Tvib, Tref)
#
#    return s0.compare_with(s_resc, spectra_only=True, plot=plot)


if __name__ == "__main__":

    from radis.test.lbl.test_overp import run_testcases

    print("test_overp.py:", run_testcases(plot=True))
