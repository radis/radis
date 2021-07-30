import pytest

from radis.lbl.factory import SpectrumFactory
from radis.test.utils import setup_test_line_databases


@pytest.mark.fast
def test_print_perf_profile(verbose=True, *args, **kwargs):

    setup_test_line_databases()  # add HITRAN-CO-TEST in ~/radis.json if not there

    sf = SpectrumFactory(
        wavelength_min=4400,
        wavelength_max=4800,
        mole_fraction=0.01,
        cutoff=1e-25,
        wstep="auto",
        isotope=[1],
        db_use_cached=True,
        self_absorption=True,
        verbose=verbose,
    )

    sf.warnings["MissingSelfBroadeningWarning"] = "ignore"
    sf.warnings["NegativeEnergiesWarning"] = "ignore"
    sf.warnings["HighTemperatureWarning"] = "ignore"

    sf.load_databank("HITRAN-CO-TEST")

    s1 = sf.eq_spectrum(300, pressure=1)

    # printing performance profiler from Spectrum object
    s1.print_perf_profile()

    # printing performance profiler from SpectrumFactory object
    sf.print_perf_profile()
