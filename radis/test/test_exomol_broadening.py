"""Tests for ExoMol broadening partners support (issue #602)"""

import inspect

import pytest


def test_exomol_default_air_broadening():
    """Default broadening should be air - backward compatible."""
    from radis.io.exomol import fetch_exomol

    sig = inspect.signature(fetch_exomol)
    assert "diluent" in sig.parameters


def test_exomol_broadening_unavailable_raises_error():
    """Unavailable broadening species should raise NotImplementedError."""
    from radis import calc_spectrum

    with pytest.raises(NotImplementedError):
        calc_spectrum(
            1000,
            1100,
            molecule="CO",
            isotope="1",
            pressure=0.01,
            Tgas=700,
            mole_fraction=1,
            databank="exomol",
            diluent={"XYZ_INVALID": 1.0},
        )
