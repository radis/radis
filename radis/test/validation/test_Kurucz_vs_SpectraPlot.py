# -*- coding: utf-8 -*-
"""
Comparing an atomic spectrum for the Kurucz database with that generated from SpectraPlot (https://spectraplot.com/)
1. `test_Kurucz_vs_NISTandSpectraplot` uses `SpectrumFactory` directly, works in nm, tests the absorption coefficient, and uses Pandas
2. `test_Kurucz_vs_NISTandSpectraplot_4000` uses `calc_spectrum`, with values for most parameters largely equivalent to the defaults in SpectraPlot, works in cm-1, tests the integral under the emission graphs, uses Vaex, and produces verbose output
Note that SpectraPlot appears to use the 'Observed Wavelength Air' given by NIST, whereas the wavelengths in the Kurucz database seem to equate to the more precise 'Ritz Wavelength Air' in NIST, so there may be a relative shift between the spectrum from RADIS and that from SpectraPlot
"""
import numpy as np
import pytest

import radis
from radis import Spectrum, SpectrumFactory, calc_spectrum, plot_diff
from radis.misc.utils import NotInstalled, not_installed_vaex_args
from radis.test.utils import getValidationCase

try:
    import vaex
except ImportError:
    vaex = NotInstalled(*not_installed_vaex_args)


def func1(**kwargs):
    """An example implementing the default broadening formula and values of SpectraPlot"""
    # print(kwargs.keys())
    # print(kwargs['df'].columns)
    return 0.1 * (296 / kwargs["Tgas"]) ** 0.8, None


@pytest.mark.needs_connection
def test_Kurucz_vs_NISTandSpectraplot(plot=True, verbose=True):
    def broad_arbitrary(**kwargs):
        """An arbitrary broadening formula in SpectraPlot (https://spectraplot.com/)"""
        return 1 * (296 / kwargs["Tgas"]) ** 0.8, None

    #%% Employ the same inputs than in example file 'spectraplot_O_10000K.txt'
    sf = SpectrumFactory(
        wavelength_min=777,
        wavelength_max=778,
        wstep=0.001,
        species="O_I",
        optimization="simple",
        path_length=1,  # cm
        pressure=1,  # atm
        verbose=0,
        lbfunc=broad_arbitrary,
    )
    sf.fetch_databank("kurucz", parfuncfmt="kurucz")
    # sf.load_databank('Kurucz-O_I', drop_columns=[], load_columns='all')
    s_RADIS = sf.eq_spectrum(Tgas=10000, name="Kurucz by RADIS")
    # s_RADIS.plot("radiance_noslit", wunit="cm-1")

    #%% Experimental spectrum
    L = 1  # cm - Input in SpectraPlot software
    raw_data = np.loadtxt(
        getValidationCase("spectraplot_O_10000K.txt"), delimiter=",", skiprows=1
    )
    s_SpectraPlot = Spectrum(
        {
            "wavenumber": raw_data[:, 0],
            "abscoeff": raw_data[:, 1]
            / L,  # spectraplot outputs absorbance; abscoef = absorbance / L
        },
        units={"abscoeff": "cm-1"},
        wunit="cm-1",
        name="NIST by SpectraPlot",
    )
    # s_SpectraPlot.plot('abscoeff', wunit='cm-1')

    if plot:
        plot_diff(s_RADIS, s_SpectraPlot, "abscoeff", wunit="nm")
    A_RADIS = s_RADIS.get_integral("abscoeff", wunit="nm")
    A_SpectraPlot = s_SpectraPlot.get_integral("abscoeff", wunit="nm")

    if verbose:
        print(
            f"Ratio of area under abscoef ('k') is A_RADIS/A_SpectraPlot = {A_RADIS/A_SpectraPlot:.3f}"
        )

    assert np.isclose(A_RADIS, A_SpectraPlot, rtol=1e-2)


@pytest.mark.needs_connection
@pytest.mark.skipif(isinstance(vaex, NotInstalled), reason="Vaex not available")
def test_Kurucz_vs_NISTandSpectraplot_4000(plot=True, verbose=True):
    initial_engine = radis.config["DATAFRAME_ENGINE"]
    radis.config["DATAFRAME_ENGINE"] = "vaex"

    w, I = np.loadtxt(
        getValidationCase("spectraplot_O_4000K.csv"),
        skiprows=1,
        delimiter=",",
        unpack=True,
    )
    I = np.where(I == 0, 1e-99, I)

    s_SpectraPlot = Spectrum.from_array(
        w, I, quantity="radiance_noslit", wunit="cm-1", Iunit="μW/cm2/sr/cm-1"
    )

    s_RADIS = calc_spectrum(
        12850,
        12870,  # from this up to 13120 is largely 0
        species="O",  # should be converted to O_I by radis.db.classes.to_conventional_name
        Tgas=4000,  # K
        databank="kurucz",
        pressure=1.01325,
        path_length=15,
        lbfunc=func1,
        warnings={"AccuracyError": "ignore", "AccuracyWarning": "ignore"},
        verbose=2,
    )

    # s_RADIS.plot("radiance_noslit", wunit="cm-1")
    if plot:
        plot_diff(
            s_SpectraPlot, s_RADIS, label1="SpectraPlot", label2="Kurucz"
        )  # , method='ratio')

    I_RADIS = s_RADIS.get_integral(
        "radiance_noslit", wunit="cm-1", Iunit="mW/cm2/sr/cm-1"
    )  # , return_units=True)
    I_SpectraPlot = s_SpectraPlot.get_integral(
        "radiance_noslit", wunit="cm-1", Iunit="mW/cm2/sr/cm-1"
    )  # , return_units=True)
    # print(I_RADIS, I_SpectraPlot, I_RADIS - I_SpectraPlot, (I_RADIS - I_SpectraPlot)/I_SpectraPlot)
    if verbose:
        print(
            f"Ratio of area under emission (radiance) is I_RADIS/I_SpectraPlot = {I_RADIS/I_SpectraPlot}"
        )

    assert np.isclose(I_RADIS, I_SpectraPlot, rtol=1.4e-2)

    radis.config["DATAFRAME_ENGINE"] = initial_engine


if __name__ == "__main__":
    test_Kurucz_vs_NISTandSpectraplot()
    test_Kurucz_vs_NISTandSpectraplot_4000()
