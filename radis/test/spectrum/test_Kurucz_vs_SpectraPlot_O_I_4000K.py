from radis.spectrum.spectrum import Spectrum
from radis import calc_spectrum, plot_diff
import numpy as np
from radis.test.utils import getTestFile

def func1(**kwargs):
    """An example implementing the default broadening formula and values of SpectraPlot (https://spectraplot.com/)"""
    # print(kwargs.keys())
    # print(kwargs['df'].columns)
    return 0.1*(296/kwargs['Tgas'])**0.8, None

def Kurucz_vs_NISTandSpectraplot(verbose=True):
    w, I = np.loadtxt(getTestFile("spectraplot_O_4000K.csv"), skiprows=1, delimiter=',', unpack=True)
    I = np.where(I==0, 1e-99, I)

    s_SpectraPlot = Spectrum.from_array(w, I, quantity='radiance_noslit', wunit='cm-1', Iunit='μW/cm2/sr/cm-1')

    s_RADIS = calc_spectrum(
        12850,
        12870,# from this up to 13120 is largely 0
        species="O_I",  # Enter species name
        Tgas=4000,  # K
        databank="kurucz",
        pressure=1.01325,
        path_length=15,
        lbfunc=func1,
        warnings={"AccuracyError": "ignore", "AccuracyWarning": "ignore"},
    )

    #s_RADIS.plot("radiance_noslit", wunit="cm-1")
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
    Kurucz_vs_NISTandSpectraplot()