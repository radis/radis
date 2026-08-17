"""
Test automatic linestrength cutoff via cutoff_error parameter.
Addresses GitHub issue #268.
"""
import pytest

from radis import SpectrumFactory


def test_cutoff_error_pandas():
    """Test that cutoff_error automatically sets cutoff for pandas dataframes."""
    sf = SpectrumFactory(
        wavenum_min=2000,
        wavenum_max=2300,
        molecule="CO",
        isotope="1",
        pressure=0.01,
        mole_fraction=1,
        path_length=1,
        cutoff=1e-27,
        cutoff_error=2,
        wstep="auto",
        verbose=0,
        warnings={"AccuracyError": "ignore", "LinestrengthCutoffWarning": "ignore"},
    )
    sf.fetch_databank("hitran")
    s = sf.eq_spectrum(Tgas=300)
    assert s is not None
    # cutoff should have been auto-adjusted from 1e-27
    assert sf.params.cutoff != 1e-27
    assert sf.params.cutoff_error == 2


def test_cutoff_error_invalid():
    """Test that invalid cutoff_error raises ValueError."""
    sf = SpectrumFactory(
        wavenum_min=2000,
        wavenum_max=2300,
        molecule="CO",
        isotope="1",
        pressure=0.01,
        mole_fraction=1,
        path_length=1,
        cutoff=1e-27,
        cutoff_error=150,  # invalid: > 100
        wstep="auto",
        verbose=0,
        warnings={"AccuracyError": "ignore"},
    )
    sf.fetch_databank("hitran")
    with pytest.raises(ValueError, match="cutoff_error must be between"):
        sf.eq_spectrum(Tgas=300)


def test_cutoff_error_vaex():
    """Test that cutoff_error works with Vaex dataframes using log-spaced binning."""
    pytest.importorskip("vaex")
    sf = SpectrumFactory(
        wavenum_min=2000,
        wavenum_max=2300,
        molecule="CO",
        isotope="1",
        pressure=0.01,
        mole_fraction=1,
        path_length=1,
        cutoff=1e-27,
        cutoff_error=2,
        wstep="auto",
        verbose=0,
        warnings={"AccuracyError": "ignore", "LinestrengthCutoffWarning": "ignore"},
    )
    sf.fetch_databank("hitran", load_columns="vaex")
    s = sf.eq_spectrum(Tgas=300)
    assert s is not None
    assert sf.params.cutoff != 1e-27
    assert sf.params.cutoff_error == 2


if __name__ == "__main__":
    test_cutoff_error_pandas()
    print("All tests passed!")
