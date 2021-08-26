# -*- coding: utf-8 -*-
"""
Benchmark RADIS performances against a HAPI equilibrium simulation of Methane

Notes
-----

We try to use the same Computation parameters when possible, in particular
the :py:attr:`~radis.lbl.factory.SpectrumFactory.broadening_max_width` of RADIS
is defined in terms of WavenumberWingHW/HWHM in [HAPI]_



References
----------

[HAPI]_ article, Table 7, case Methane III with SCPF, 50 HWHMs  :

    >>> Calculation time in article: 317 s


Typical results on an XPS 15 laptop here::

    >>> Calculated with HAPI in 157.41s
    >>> Calculated with RADIS in 1.65s


-------------------------------------------------------------------------------

"""

from os.path import dirname, join
from time import time

from hapi import ISO_ID, absorptionCoefficient_Voigt, db_begin, fetch_by_ids, tableList

from radis import Spectrum, SpectrumFactory, plot_diff

if __name__ == "__main__":

    benchmark_line_brd_ratio = 50  # “WavenumberWingHW”/HWHMs
    dnu = 0.01  # step in HAPI Benchmark article
    molecule = "CH4"
    wavenum_min = 0.001
    wavenum_max = 11505
    pressure_bar = 1.01315
    T = 296
    isotopes = [1, 2, 3, 4]

    sf = SpectrumFactory(
        wavenum_min=wavenum_min,
        wavenum_max=wavenum_max,
        isotope=isotopes,  #'all',
        verbose=2,
        wstep=dnu,  # depends on HAPI benchmark.
        cutoff=1e-23,
        truncation=2.865,  # Corresponds to WavenumberWingHW/HWHM=50 in HAPI
        molecule=molecule,
        optimization=None,
    )
    sf.fetch_databank("hitran")

    s = sf.eq_spectrum(Tgas=T, pressure=pressure_bar)
    s.name = "RADIS ({0:.1f}s)".format(s.conditions["calculation_time"])

    # Print our HWHM for comparison (a posteriori)
    print(("HWHM max {0:.2f} cm-1".format(sf.df1.hwhm_voigt.max())))
    print(
        (
            "WavenumberWingHW/HWHM",
            int(sf.params.truncation / (sf.df1.hwhm_voigt.max())),
        )
    )
    assert (
        int(sf.params.truncation / (sf.df1.hwhm_voigt.max()))
    ) == benchmark_line_brd_ratio

    # %% Run HAPI
    print("Now calculating with HAPI")

    # Generate HAPI database locally
    db_begin(join(dirname(__file__), __file__.replace(".py", "_HAPIdata")))
    if not molecule in tableList():  # only if data not downloaded already
        mol_iso_ids = [k for k in ISO_ID if ISO_ID[k][-1] == molecule]
        mol_iso_ids = [k for k in mol_iso_ids if ISO_ID[k][1] in isotopes]
        fetch_by_ids(molecule, mol_iso_ids, wavenum_min, wavenum_max)

    # Calculate with HAPI
    def calc_hapi():
        nu, coef = absorptionCoefficient_Voigt(
            SourceTables=molecule,
            Environment={
                "T": T,
                "p": pressure_bar / 1.01315,
            },  # K  # atm
            GammaL="gamma_self",
            WavenumberStep=dnu,
            HITRAN_units=False,
        )
        return nu, coef

    t0 = time()
    nu, coef = calc_hapi()
    t0 = time() - t0
    print(("Calculated with HAPI in {0:.2f}s".format(t0)))

    s_hapi = Spectrum.from_array(nu, coef, "abscoeff", wunit="cm-1", unit="cm-1")
    s_hapi.name = "HAPI ({0:.1f}s)".format(t0)

    plot_diff(s_hapi, s, "abscoeff")

    print(
        ("Calculated with RADIS in {0:.2f}s".format(s.conditions["calculation_time"]))
    )
    print(("Number of lines in RADIS:", len(sf.df0)))

#    plt.savefig('radis_vs_hapi_test_large_ch4.pdf')
