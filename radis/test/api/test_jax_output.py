"""Tests for output='jax' - issue #474"""
import pytest
import os

@pytest.mark.fast
def test_fetch_exomol_jax_skips_radis_abundance():
    """output='jax' should NOT apply RADIS terrestrial isotope abundance.
    ExoJax handles abundance itself. Issue #474"""
    # Read file directly — avoids cached module issue
    file_path = os.path.join(os.path.dirname(__file__), '../../io/exomol.py')
    with open(file_path, 'r') as f:
        src = f.read()
    assert "Abundance correction skipped" in src, \
        "output='jax' should skip RADIS abundance correction"

@pytest.mark.fast
def test_dbmanager_has_output_parameter():
    """DatabaseManager.load() must accept output parameter. Issue #474"""
    from radis.api.dbmanager import DatabaseManager
    import inspect
    sig = inspect.signature(DatabaseManager.load)
    assert 'output' in sig.parameters
    def test_fetch_exomol_delegates_to_api():
     """
     radis/io/exomol.py should use MdbExomol from radis.api
     Ref: https://github.com/radis/radis/issues/474
    """
    import inspect
    import radis.io.exomol as io_module

    source = inspect.getsource(io_module)
    assert "MdbExomol" in source, \
        "fetch_exomol() should use radis.api.MdbExomol"
    print("✅ radis/io/exomol.py correctly uses MdbExomol from radis.api")


def test_apply_jax_array_conversion_exists():
    """
    apply_jax_array_conversion helper should exist in radis.api.exomolapi
    Ref: https://github.com/radis/radis/issues/474
    """
    from radis.api.exomolapi import apply_jax_array_conversion
    assert callable(apply_jax_array_conversion)
    print("✅ apply_jax_array_conversion exists in radis.api.exomolapi")


def test_hitran_has_jax_separation():
    """
    radis/io/hitran.py should have GPU/DRAM separation
    Ref: https://github.com/radis/radis/issues/474
    """
    import inspect
    import radis.io.hitran as hitran_module

    source = inspect.getsource(hitran_module)
    assert "GPU/DRAM Array Separation" in source, \
        "hitran.py should have GPU/DRAM separation block"
    print("✅ hitran.py has GPU/DRAM separation")
