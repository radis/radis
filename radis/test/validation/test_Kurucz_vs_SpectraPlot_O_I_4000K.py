# -*- coding: utf-8 -*-
"""
Comparing an atomic spectrum for the Kurucz database with that generated from SpectraPlot (https://spectraplot.com/), specifically:
using calc_spectrum, with values for the parameters largely equivalent to the defaults in SpectraPlot, working in cm-1, and testing the integral under the emission graphs, and using vaex rather than pandas, and verbose output
"""

from radis.spectrum.spectrum import Spectrum
from radis import calc_spectrum, plot_diff
import numpy as np
from radis.test.utils import getValidationCase
import pytest
import radis

radis.config['DATAFRAME_ENGINE'] = 'vaex'

def func1(**kwargs):
    """An example implementing the default broadening formula and values of SpectraPlot"""
    # print(kwargs.keys())
    # print(kwargs['df'].columns)
    return 0.1*(296/kwargs['Tgas'])**0.8, None

@pytest.mark.needs_connection
def test_Kurucz_vs_NISTandSpectraplot_4000(plot=True, verbose=True):
    w, I = np.loadtxt(getValidationCase("spectraplot_O_4000K.csv"), skiprows=1, delimiter=',', unpack=True)
    I = np.where(I==0, 1e-99, I)

    s_SpectraPlot = Spectrum.from_array(w, I, quantity='radiance_noslit', wunit='cm-1', Iunit='μW/cm2/sr/cm-1')

    s_RADIS = calc_spectrum(
        12850,
        12870,# from this up to 13120 is largely 0
        species="O",  # should be converted to O_I by radis.db.classes.to_conventional_name
        Tgas=4000,  # K
        databank="kurucz",
        pressure=1.01325,
        path_length=15,
        lbfunc=func1,
        warnings={"AccuracyError": "ignore", "AccuracyWarning": "ignore"},
        verbose=2
    )

    #s_RADIS.plot("radiance_noslit", wunit="cm-1")
    if plot:
        plot_diff(s_SpectraPlot, s_RADIS, label1='SpectraPlot', label2='Kurucz')#, method='ratio')

    I_RADIS = s_RADIS.get_integral("radiance_noslit", wunit="cm-1", Iunit='mW/cm2/sr/cm-1')#, return_units=True)
    I_SpectraPlot = s_SpectraPlot.get_integral("radiance_noslit", wunit="cm-1", Iunit='mW/cm2/sr/cm-1')#, return_units=True)
    #print(I_RADIS, I_SpectraPlot, I_RADIS - I_SpectraPlot, (I_RADIS - I_SpectraPlot)/I_SpectraPlot)
    if verbose:
        print(
            f"Ratio of area under emission (radiance) is I_RADIS/I_SpectraPlot = {I_RADIS/I_SpectraPlot}"
        )

    assert np.isclose(I_RADIS, I_SpectraPlot, rtol=1.4e-2)

if __name__ == "__main__":
    test_Kurucz_vs_NISTandSpectraplot_4000()