# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 14:29:42 2017

@author: erwan

Summary
-------

Reproducing the validation case of Klarenaar 2017 [1]_, who calculated a transmittance
spectrum from the initial data of Dang 1982 [2]_, with a 1 rotational temperature +
3 vibrational temperature (Treanor distributions) model

CO2 Energies are calculated from Dunham developments in an uncoupled harmonic
oscillator - rigid rotor model

References
----------

.. [1] Klarenaar et al 2017, "Time evolution of vibrational temperatures in a CO2 glow
       discharge measured with infrared absorption spectroscopy" doi/10.1088/1361-6595/aa902e

.. [2] Dang et al 1982, "Detailed vibrational population distributions in a CO2 laser
        discharge as measured with a tunable diode laser" doi/10.1007/BF00694640


-------------------------------------------------------------------------------


"""

from os.path import join

from radis import SpectrumFactory
from radis.misc.printer import printm
from radis.spectrum import Spectrum, get_residual, plot_diff
from radis.test.utils import getValidationCase, setup_test_line_databases


def test_klarenaar_validation_case(
    verbose=True, plot=False, warnings=True, *args, **kwargs
):
    """Reproduce the Klarenaar 2018 validation case, as given in the
    [RADIS-2018]_ article.

    References
    ----------

    Klarenaar et al, "Time evolution of vibrational temperatures in a CO 2 glow
    discharge measured with infrared absorption spectroscopy", doi 10.1088/1361-6595/aa902e,
    and the references there in.

    """

    setup_test_line_databases()

    # %% Data from Dang, adapted by Klarenaar
    s_exp = Spectrum.from_txt(
        getValidationCase(
            join(
                "test_CO2_3Tvib_vs_klarenaar_data", "klarenaar_2017_digitized_data.csv"
            )
        ),
        "transmittance_noslit",
        wunit="cm-1",
        unit="",
        delimiter=",",
        name="Klarenaar 2017",
    )

    # %% Calculate Klarenaar test case conditions

    sf = SpectrumFactory(
        2284.2,
        2284.6,
        wstep=0.001,  # cm-1
        pressure=20 * 1e-3,  # bar
        cutoff=1e-25,
        molecule="CO2",
        isotope="1,2",
        path_length=10,  # cm-1
        # warning! 10% in mass fraction -> less in mole fraction
        mole_fraction=0.1 * 28.97 / 44.07,
        truncation=0.5,  # cm-1
        medium="vacuum",
        export_populations="vib",
    )
    sf.warnings["MissingSelfBroadeningWarning"] = "ignore"
    #        sf.load_databank('HITEMP-CO2-DUNHAM')
    sf.load_databank("HITEMP-CO2-TEST")

    # Calculate with Klarenaar fitted values
    T12 = 517
    T3 = 2641
    Trot = 491

    s = sf.non_eq_spectrum(
        (T12, T12, T3), Trot, Ttrans=Trot, vib_distribution="treanor", name="RADIS"
    )

    if plot:
        plot_diff(s, s_exp, "transmittance_noslit")
    #        plt.savefig('test_CO2_3Tvib_vs_klarenaar.png')

    assert get_residual(s, s_exp, "transmittance_noslit", ignore_nan=True) < 0.003

    return True


if __name__ == "__main__":
    printm(
        "test_CO2_3Tvib_vs_klarenaar:",
        test_klarenaar_validation_case(verbose=True, plot=True),
    )
