# -*- coding: utf-8 -*-
"""
Example using GPU sliders and GPU calculation with :py:meth:`~radis.lbl.SpectrumFactory.eq_spectrum_gpu`

Use ``emulate=True`` to run the GPU code on CPU, and ``emulate=False`` (default)
to run it with the full power of the GPU

"""

from radis import SpectrumFactory
from radis.test.utils import getTestFile
from radis.tools.database import load_spec
from radis.tools.plot_tools import ParamRange

my_file = getTestFile("CO2_measured_spectrum_4-5um.spec")  # for the example here
s_exp = load_spec(my_file)

s_exp.crop(4120, 4790).plot(Iunit="mW/cm2/sr/nm")


sf = SpectrumFactory(
    2100,
    2400,  # cm-1
    molecule="CO2",
    isotope="1,2,3",
    wstep=0.002,
)

sf.fetch_databank("hitemp")

s = sf.eq_spectrum_gpu_explore(
    var="radiance",
    Tgas=ParamRange(300.0, 2500.0, 1100.0),  # K
    pressure=ParamRange(0.1, 2, 1),  # bar
    mole_fraction=ParamRange(0, 1, 0.8),
    path_length=ParamRange(0, 1, 0.2),  # cm
    slit_FWHM=ParamRange(0, 1.5, 0.24),  # cm-1
    emulate=True,  # runs on CPU
    plotkwargs={"nfig": "same", "wunit": "nm"},
)
print(s)
