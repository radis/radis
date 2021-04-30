# -*- coding: utf-8 -*-
"""
Scripts used to generate some test files for RADIS


-------------------------------------------------------------------------------

"""


from radis.test.utils import setup_test_line_databases

if __name__ == "__main__":

    # %% Generate carbon monoxide files

    from radis import SpectrumFactory

    setup_test_line_databases()  # add HITRAN-CO-TEST in ~/radis.json if not there

    Tgas = 1500
    sf = SpectrumFactory(
        wavelength_min=4400,
        wavelength_max=4800,
        mole_fraction=0.01,
        #                         path_length=0.1,
        cutoff=1e-25,
        wstep=0.005,
        isotope=[1],
        db_use_cached=True,
        self_absorption=True,
        verbose=False,
    )
    sf.load_databank("HITRAN-CO-TEST")
    s1 = sf.non_eq_spectrum(Tgas, Tgas, path_length=0.01)
    s1.store("CO_Tgas1500K_mole_fraction0.01.spec", compress=True)

    s2 = sf.non_eq_spectrum(Tgas, Tgas, path_length=0.01, mole_fraction=0.5)
    s2.store("CO_Tgas1500K_mole_fraction0.5.spec", compress=True)
