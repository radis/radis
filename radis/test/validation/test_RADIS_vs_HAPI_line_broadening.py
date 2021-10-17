# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 21:10:31 2017

@author: erwan

Validation case
------------

Compare line broadening calculated in RADIS vs calculated in HAPI, using the
same database (HITRAN 2016)

"""

import shutil
from os.path import dirname, exists, join

import pytest

from radis import Spectrum, SpectrumFactory
from radis.db.classes import get_molecule_identifier
from radis.misc.printer import printm
from radis.phys.convert import nm2cm
from radis.test.utils import setup_test_line_databases


@pytest.mark.fast
@pytest.mark.needs_connection  # ignored by pytest with argument -m "not needs_http_connection"
def test_line_broadening(rtol=1e-3, verbose=True, plot=False, *args, **kwargs):
    r"""
    Plot absorption coefficient (cm-1) of CO at high temperature (2000 K) with
    RADIS, and compare with calculations from HAPI using the HITRAN database

    Notes
    -----

    In this example no data is needed. Everything is downloaded from the HITRAN
    database directly using either the HAPI ``fetch`` function, or the RADIS
    :meth:`~neq.spec.factory.fetch_databank` method.

    """

    from hapi import (
        absorptionCoefficient_Voigt,
        db_begin,
        fetch,
        tableList,
        transmittanceSpectrum,
    )

    setup_test_line_databases()  # add HITRAN-CO-TEST in ~/radis.json if not there

    # Conditions
    molecule = "CO2"
    mol_id = get_molecule_identifier(molecule)
    iso = 1
    T = 1500
    p = 0.1
    L = 0.1
    #    M = 0.001           # mole fraction  (dont know where to put that )
    dnu = 0.0001
    wmin = nm2cm(4372.69 + 0.2)  # cm-1
    wmax = nm2cm(4372.69 - 0.2)  # cm-1
    truncation = 0.25  # cm-1
    neighbour_lines = 0.5  # cm-1

    # %% HITRAN calculation
    # -----------

    # Generate HAPI database locally

    HAPIdb = join(dirname(__file__), __file__.replace(".py", "_HAPIdata"))

    def calc_hapi():
        """Calc spectrum under HAPI"""

        clean_after_run = not exists(HAPIdb) and False

        try:
            db_begin(HAPIdb)
            if not molecule in tableList():  # only if data not downloaded already
                fetch(
                    molecule,
                    mol_id,
                    iso,
                    wmin - neighbour_lines,
                    wmax + neighbour_lines,
                )
                # HAPI doesnt correct for side effects

            # Calculate with HAPI
            nu, coef = absorptionCoefficient_Voigt(
                SourceTables="CO2",
                Environment={
                    "T": T,
                    "p": p / 1.01325,
                },  # K  # atm
                WavenumberStep=dnu,
                HITRAN_units=False,
                GammaL="gamma_self",
            )
            nu, trans = transmittanceSpectrum(
                nu,
                coef,
                Environment={
                    "l": L,
                },
            )  # cm
            s_hapi = Spectrum.from_array(
                nu,
                trans,
                "transmittance_noslit",
                "cm-1",
                "1",
                conditions={"Tgas": T},
                name="HAPI",
            )

        except:
            raise

        finally:
            if clean_after_run:
                shutil.rmtree(HAPIdb)
        return s_hapi

    s_hapi = calc_hapi()

    def calc_radis():

        # %% Calculate with RADIS
        # ----------
        pl = SpectrumFactory(
            wavenum_min=wmin,
            wavenum_max=wmax,
            mole_fraction=1,
            path_length=L,
            wstep=dnu,
            molecule=molecule,
            pressure=p,
            truncation=truncation,
            cutoff=1e-23,
            isotope=iso,
        )
        pl.warnings["MissingSelfBroadeningWarning"] = "ignore"
        pl.warnings["HighTemperatureWarning"] = "ignore"
        pl.fetch_databank(
            source="hitran",
            db_use_cached=True,
        )

        s = pl.eq_spectrum(Tgas=T)  # , Ttrans=300)
        s.name = "RADIS"

        if plot:
            pl.plot_broadening()

        return s

    s = calc_radis()

    # %% Compare
    # also shrinks HAPI range to the valid one
    s_hapi.resample(s.get_wavenumber(), unit="cm-1", energy_threshold=0.1)

    save = False  # just used in the article
    if plot or save:
        from radis import plot_diff

        #        title = '{0} bar, {1} K, {2} cm'.format(p, T, L)  if save else None
        fig, [ax0, ax1] = plot_diff(
            s, s_hapi, var="transmittance_noslit", method="ratio", show=plot
        )

        ax0.annotate(
            r"[P64](00$^\mathregular{0}$0)$\rightarrow $(00$^\mathregular{0}$1)",
            (2286.945, 0.76),
            (2286.94, 0.8),
            arrowprops=dict(arrowstyle="->", facecolor="black"),
        )
        ax0.annotate(
            r"[P53](01$^\mathregular{1}$0)$\rightarrow $(01$^\mathregular{1}$1)",
            (2286.9, 0.78),
            (2286.9, 0.82),
            arrowprops=dict(arrowstyle="->", facecolor="black"),
            horizontalalignment="right",
        )
        ax1.set_ylim(0.95, 1.05)

        if save:
            fig.savefig("out/test_RADIS_vs_HAPI_line_broadening.pdf")

    # Compare integrals
    diff = abs(
        s.get_integral("transmittance_noslit")
        / s_hapi.get_integral("transmittance_noslit")
        - 1
    )
    b = diff < rtol

    if verbose:
        printm(
            "Integral difference ({0:.2f}%) < {1:.2f}%: {2}".format(
                diff * 100, rtol * 100, b
            )
        )

    return b


if __name__ == "__main__":
    printm("test RADIS_vs_HAPI_line_broadening: ", test_line_broadening(plot=True))
