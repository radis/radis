import numpy as np
import pandas as pd


class MockSpectrumFactory:
    """Mock SpectrumFactory for unit-testing uncertainty propagation.

    Provides a minimal ``eq_spectrum()`` that returns a deterministic
    spectrum influenced by ``Tgas``.  Also exposes a ``df0`` DataFrame
    with HITRAN-like columns (``int``, ``airbrd``, ``selbrd``,
    ``Tdpair``, ``Pshft``, ``ierr``) so that line-parameter
    perturbation can be tested without a real databank.
    """

    def __init__(self, n_lines=20):
        self.input = type(
            "obj",
            (object,),
            {
                "isotope": "1",
                "mole_fraction": 0.1,
                "pressure": 1.0,
                "path_length": 10.0,
            },
        )
        self.params = type(
            "obj",
            (object,),
            {"wavenum_min": 2000, "wavenum_max": 2300},
        )
        self.conditions = {}
        self.units = {}
        self.wunit = "cm-1"

        # Build a minimal df0 with HITRAN-like columns + ierr
        rng = np.random.default_rng(0)
        self.df0 = pd.DataFrame(
            {
                "wav": np.linspace(2000, 2300, n_lines),
                "int": rng.uniform(1e-22, 1e-20, n_lines),
                "airbrd": rng.uniform(0.05, 0.10, n_lines),
                "selbrd": rng.uniform(0.08, 0.12, n_lines),
                "Tdpair": rng.uniform(0.5, 0.8, n_lines),
                "Pshft": rng.uniform(-0.005, 0.005, n_lines),
                # ierr: 6-digit code — use 555500 (good codes for the
                # first 4 params, unreported for last 2)
                "ierr": np.full(n_lines, 555500, dtype=int),
            }
        )

    def eq_spectrum(self, **kwargs):
        from radis.spectrum.spectrum import Spectrum

        w = np.linspace(2000, 2300, 100)
        I = np.sin(w / 100) + kwargs.get("Tgas", 0) / 1000
        return Spectrum(
            {"radiance_noslit": (w, I)},
            wunit="cm-1",
            units={"radiance_noslit": "W/cm2/sr/cm-1"},
            conditions=kwargs,
        )
