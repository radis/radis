# -*- coding: utf-8 -*-
"""
Tests for the Spectral Uncertainty Quantification module.

Tests cover:
- HITRAN uncertainty code decoding
- GEISA uncertainty decoding
- UncertaintyModel creation and sampling
- ConfidenceBandResult statistics
- Spectrum confidence band storage
- SensitivityAnalyzer parameter registration
- get_auto_drop_columns helper
"""

import numpy as np
import pytest

from radis.tools.uncertainty import (
    HITRAN_IERR_PARAMS,
    ConfidenceBandResult,
    SensitivityAnalyzer,
    UncertaintyModel,
    UncertaintyPropagator,
    decode_geisa_uncertainty,
    decode_hitran_uncertainty,
    parse_hitran_ierr,
)

# ============================================================================
# HITRAN Uncertainty Decoder Tests
# ============================================================================


@pytest.mark.fast
def test_decode_hitran_uncertainty_all_codes():
    """Test that all 10 codes (0-9) return valid results for
    intensity (relative uncertainty)."""
    for code in range(10):
        result = decode_hitran_uncertainty(code, "intensity")
        if code == 0:
            assert result is None, f"Code 0 should be None, got {result}"
        else:
            assert isinstance(
                result, tuple
            ), f"Code {code} should return tuple, got {type(result)}"
            assert len(result) == 2
            min_val, max_val = result
            assert isinstance(min_val, float)
            # max_val can be None (for code 1) or float
            if max_val is not None:
                assert max_val > min_val


@pytest.mark.fast
def test_decode_hitran_uncertainty_position():
    """Test absolute position uncertainty decoding."""
    # Code 4 for position: 0.001 to 0.01 cm⁻¹
    result = decode_hitran_uncertainty(4, "position")
    assert result == (0.001, 0.01)

    # Code 0: unreported
    result = decode_hitran_uncertainty(0, "position")
    assert result is None


@pytest.mark.fast
def test_decode_hitran_uncertainty_relative():
    """Test relative uncertainty for intensity."""
    # Code 5: 1% to 2%
    result = decode_hitran_uncertainty(5, "intensity")
    assert result == (1.0, 2.0)

    # Code 3: 5% to 10%
    result = decode_hitran_uncertainty(3, "air_broadened")
    assert result == (5.0, 10.0)


@pytest.mark.fast
def test_decode_hitran_uncertainty_pressure_shift():
    """Test pressure shift (absolute) uncertainty."""
    result = decode_hitran_uncertainty(3, "pressure_shift")
    assert result is not None
    assert isinstance(result, tuple)


@pytest.mark.fast
def test_decode_hitran_uncertainty_invalid_code():
    """Test that invalid codes raise ValueError."""
    with pytest.raises(ValueError, match="must be 0–9"):
        decode_hitran_uncertainty(10, "intensity")
    with pytest.raises(ValueError, match="must be 0–9"):
        decode_hitran_uncertainty(-1, "intensity")


@pytest.mark.fast
def test_decode_hitran_uncertainty_invalid_parameter():
    """Test that unknown parameter names raise ValueError."""
    with pytest.raises(ValueError, match="Unknown parameter"):
        decode_hitran_uncertainty(5, "nonexistent_param")


@pytest.mark.fast
def test_decode_hitran_uncertainty_all_params():
    """Test that all valid parameter names work."""
    for param in HITRAN_IERR_PARAMS.values():
        result = decode_hitran_uncertainty(5, param)
        assert result is not None, f"Code 5 for '{param}' should not be None"


@pytest.mark.fast
def test_parse_hitran_ierr():
    """Test parsing a full 6-digit ierr code."""
    # Code 554400: position=5, intensity=5, air=4, self=4,
    #              temp_dep=0, shift=0
    result = parse_hitran_ierr(554400)
    assert isinstance(result, dict)
    assert len(result) == 6
    assert "position" in result
    assert "intensity" in result
    assert result["position"] is not None  # code 5
    assert result["intensity"] is not None  # code 5
    assert result["temperature_dependence"] is None  # code 0
    assert result["pressure_shift"] is None  # code 0


@pytest.mark.fast
def test_parse_hitran_ierr_zero_padded():
    """Test that small integers get zero-padded correctly."""
    # ierr = 50000 should be parsed as '050000'
    result = parse_hitran_ierr(50000)
    assert result["position"] is None  # first digit is 0
    assert result["intensity"] is not None  # second digit is 5


# ============================================================================
# GEISA Uncertainty Decoder Tests
# ============================================================================


@pytest.mark.fast
def test_decode_geisa_uncertainty_valid():
    """Test decoding a valid GEISA uncertainty value."""
    result = decode_geisa_uncertainty(0.005, "position")
    assert result == 0.005


@pytest.mark.fast
def test_decode_geisa_uncertainty_zero():
    """Test that zero uncertainty returns None (unreported)."""
    result = decode_geisa_uncertainty(0.0, "position")
    assert result is None


@pytest.mark.fast
def test_decode_geisa_uncertainty_nan():
    """Test that NaN uncertainty returns None."""
    result = decode_geisa_uncertainty(float("nan"), "position")
    assert result is None


@pytest.mark.fast
def test_decode_geisa_uncertainty_none():
    """Test that None input returns None."""
    result = decode_geisa_uncertainty(None, "position")
    assert result is None


# ============================================================================
# UncertaintyModel Tests
# ============================================================================


@pytest.mark.fast
def test_uncertainty_model_creation():
    """Test creating an UncertaintyModel."""
    model = UncertaintyModel()
    assert model.source is None
    assert not model.has_line_uncertainty
    assert model.condition_params == []


@pytest.mark.fast
def test_uncertainty_model_add_line_uncertainty():
    """Test enabling line parameter uncertainty."""
    model = UncertaintyModel()
    model.add_line_uncertainty(source="hitran")
    assert model.has_line_uncertainty
    assert model.line_uncertainty_source == "hitran"


@pytest.mark.fast
def test_uncertainty_model_add_line_uncertainty_invalid():
    """Test that invalid source raises ValueError."""
    model = UncertaintyModel()
    with pytest.raises(ValueError, match="Unknown uncertainty source"):
        model.add_line_uncertainty(source="invalid_db")


@pytest.mark.fast
def test_uncertainty_model_add_condition_normal():
    """Test adding a normal distribution condition."""
    model = UncertaintyModel()
    model.add_condition("Tgas", distribution="normal", std=50)
    assert "Tgas" in model.condition_params


@pytest.mark.fast
def test_uncertainty_model_add_condition_uniform():
    """Test adding a uniform distribution condition."""
    model = UncertaintyModel()
    model.add_condition("mole_fraction", distribution="uniform", low=0.08, high=0.12)
    assert "mole_fraction" in model.condition_params


@pytest.mark.fast
def test_uncertainty_model_add_condition_missing_params():
    """Test that missing distribution params raise ValueError."""
    model = UncertaintyModel()
    with pytest.raises(ValueError, match="requires 'std'"):
        model.add_condition("Tgas", distribution="normal")
    with pytest.raises(ValueError, match="requires 'low'"):
        model.add_condition("Tgas", distribution="uniform", low=0.5)  # missing high


@pytest.mark.fast
def test_uncertainty_model_sample_condition_normal():
    """Test sampling from a normal distribution."""
    model = UncertaintyModel()
    model.add_condition("Tgas", distribution="normal", std=50)
    samples = model.sample_condition(
        "Tgas", nominal_value=1500, n_samples=1000, seed=42
    )
    assert len(samples) == 1000
    assert abs(np.mean(samples) - 1500) < 10  # should be ~1500
    assert abs(np.std(samples) - 50) < 5  # should be ~50


@pytest.mark.fast
def test_uncertainty_model_sample_condition_uniform():
    """Test sampling from a uniform distribution."""
    model = UncertaintyModel()
    model.add_condition("mole_fraction", distribution="uniform", low=0.08, high=0.12)
    samples = model.sample_condition(
        "mole_fraction", nominal_value=0.1, n_samples=1000, seed=42
    )
    assert len(samples) == 1000
    assert np.all(samples >= 0.08)
    assert np.all(samples <= 0.12)


@pytest.mark.fast
def test_uncertainty_model_sample_condition_not_defined():
    """Test that sampling undefined param raises KeyError."""
    model = UncertaintyModel()
    with pytest.raises(KeyError, match="No uncertainty defined"):
        model.sample_condition("Tgas", nominal_value=1500)


@pytest.mark.fast
def test_uncertainty_model_invalid_distribution():
    """Test that unknown distribution type raises ValueError."""
    model = UncertaintyModel()
    with pytest.raises(ValueError, match="Unknown distribution"):
        model.add_condition("Tgas", distribution="cauchy", scale=50)


# ============================================================================
# ConfidenceBandResult Tests
# ============================================================================


@pytest.mark.fast
def test_confidence_band_result():
    """Test ConfidenceBandResult statistics."""
    rng = np.random.default_rng(42)
    wavenumber = np.linspace(2000, 2300, 100)
    # Generate 50 "spectra" — each a smooth curve + noise
    base_spectrum = np.sin(wavenumber / 100)
    spectra_array = base_spectrum + rng.normal(0, 0.1, (50, 100))

    result = ConfidenceBandResult(
        wavenumber=wavenumber,
        spectra_array=spectra_array,
        quantity="radiance_noslit",
        confidence_level=0.95,
    )

    assert result.n_samples == 50
    assert result.mean.shape == (100,)
    assert result.std.shape == (100,)
    assert result.lower.shape == (100,)
    assert result.upper.shape == (100,)

    # lower < mean < upper should hold on average
    assert np.all(result.lower <= result.upper)

    # Get confidence band
    w, lower, upper = result.get_confidence_band()
    assert np.array_equal(w, wavenumber)

    # Get stats
    stats = result.get_stats()
    assert "mean" in stats
    assert "std" in stats
    assert "median" in stats
    assert stats["n_samples"] == 50
    assert stats["confidence_level"] == 0.95


# ============================================================================
# Spectrum Confidence Band Tests
# ============================================================================


@pytest.mark.fast
def test_spectrum_confidence_band_storage():
    """Test storing and retrieving confidence bands in Spectrum."""
    from radis.spectrum.spectrum import Spectrum

    w = np.linspace(2000, 2300, 100)
    I = np.sin(w / 100)

    s = Spectrum(
        {"radiance_noslit": (w, I)},
        units={"radiance_noslit": "mW/cm2/sr/cm-1"},
        wunit="cm-1",
    )

    # Initially no confidence bands
    assert not s.has_confidence_band()
    assert not s.has_confidence_band("radiance_noslit")

    # Set confidence band
    lower = I - 0.1
    upper = I + 0.1
    s.set_confidence_band(
        "radiance_noslit",
        lower=lower,
        upper=upper,
        confidence_level=0.95,
        n_samples=200,
        std=np.full(100, 0.1),
    )

    assert s.has_confidence_band()
    assert s.has_confidence_band("radiance_noslit")
    assert not s.has_confidence_band("transmittance")

    # Get confidence band
    w_out, lower_out, upper_out = s.get_confidence_band("radiance_noslit")
    assert np.allclose(lower_out, lower)
    assert np.allclose(upper_out, upper)

    # Get stats
    stats = s.get_uncertainty_stats("radiance_noslit")
    assert stats["confidence_level"] == 0.95
    assert stats["n_samples"] == 200


@pytest.mark.fast
def test_spectrum_confidence_band_not_stored():
    """Test that accessing non-existent band raises KeyError."""
    from radis.spectrum.spectrum import Spectrum

    w = np.linspace(2000, 2300, 100)
    I = np.sin(w / 100)
    s = Spectrum(
        {"radiance_noslit": (w, I)},
        units={"radiance_noslit": "mW/cm2/sr/cm-1"},
        wunit="cm-1",
    )

    with pytest.raises(KeyError, match="No confidence band"):
        s.get_confidence_band("radiance_noslit")
    with pytest.raises(KeyError, match="No uncertainty stats"):
        s.get_uncertainty_stats("radiance_noslit")


# ============================================================================
# SensitivityAnalyzer Tests
# ============================================================================


@pytest.mark.fast
def test_sensitivity_analyzer_add_parameter():
    """Test adding parameters to SensitivityAnalyzer."""
    sa = SensitivityAnalyzer(sf=None)
    sa.add_parameter("Tgas", bounds=[1200, 1800])
    sa.add_parameter("mole_fraction", bounds=[0.05, 0.15])
    assert sa.parameter_names == ["Tgas", "mole_fraction"]


@pytest.mark.fast
def test_sensitivity_analyzer_invalid_bounds():
    """Test invalid bounds raise ValueError."""
    sa = SensitivityAnalyzer(sf=None)
    with pytest.raises(ValueError, match="low, high"):
        sa.add_parameter("Tgas", bounds=[1800])
    with pytest.raises(ValueError, match="Lower bound"):
        sa.add_parameter("Tgas", bounds=[1800, 1200])


@pytest.mark.fast
def test_sensitivity_analyzer_stub():
    """Test that sobol_analysis raises ImportError when SALib is missing."""
    import sys
    import types

    # Temporarily block SALib by replacing all SALib modules with a
    # sentinel module whose __getattr__ raises ImportError, and by
    # also inserting a None for the submodules we import.
    saved = {}
    keys_to_block = [k for k in sys.modules if k == "SALib" or k.startswith("SALib.")]
    for k in keys_to_block:
        saved[k] = sys.modules.pop(k)

    # Insert a blocker so 'from SALib.analyze import sobol' fails
    blocker = types.ModuleType("SALib")
    blocker.__path__ = []  # make it a package

    def _raise(*a, **kw):
        raise ImportError("blocked")

    blocker.__getattr__ = _raise
    sys.modules["SALib"] = blocker
    sys.modules["SALib.analyze"] = None
    sys.modules["SALib.analyze.sobol"] = None
    sys.modules["SALib.sample"] = None
    sys.modules["SALib.sample.sobol"] = None

    try:
        sa = SensitivityAnalyzer(sf=None)
        sa.add_parameter("Tgas", bounds=[1200, 1800])
        with pytest.raises(ImportError, match="SALib"):
            sa.sobol_analysis()
    finally:
        # Remove blockers
        for k in list(sys.modules.keys()):
            if k == "SALib" or k.startswith("SALib."):
                del sys.modules[k]
        # Restore original modules
        sys.modules.update(saved)


# ============================================================================
# get_auto_drop_columns Tests
# ============================================================================


@pytest.mark.fast
def test_get_auto_drop_columns_default():
    """Test default column dropping includes ierr."""
    from radis.lbl.loader import get_auto_drop_columns

    cols = get_auto_drop_columns("hitran", None)
    assert "ierr" in cols


@pytest.mark.fast
def test_get_auto_drop_columns_preserve_uncertainty():
    """Test that preserve_uncertainty=True keeps ierr."""
    from radis.lbl.loader import get_auto_drop_columns

    cols = get_auto_drop_columns("hitran", None, preserve_uncertainty=True)
    assert "ierr" not in cols


@pytest.mark.fast
def test_get_auto_drop_columns_hitemp_preserve():
    """Test preserve_uncertainty for hitemp format."""
    from radis.lbl.loader import get_auto_drop_columns

    cols_default = get_auto_drop_columns("hitemp", None)
    assert "ierr" in cols_default

    cols_preserved = get_auto_drop_columns("hitemp", None, preserve_uncertainty=True)
    assert "ierr" not in cols_preserved

    # Other columns should still be dropped
    assert "iref" in cols_preserved
    assert "lmix" in cols_preserved


@pytest.mark.fast
def test_get_auto_drop_columns_geisa():
    """Test that GEISA format has no ierr to drop."""
    from radis.lbl.loader import get_auto_drop_columns

    cols = get_auto_drop_columns("geisa", None)
    assert "ierr" not in cols

    # preserve_uncertainty should be a no-op for GEISA
    cols_preserved = get_auto_drop_columns("geisa", None, preserve_uncertainty=True)
    assert cols == cols_preserved


# ============================================================================
# UncertaintyPropagator Tests (unit-level, no actual SF)
# ============================================================================


@pytest.mark.fast
def test_uncertainty_propagator_creation():
    """Test creating an UncertaintyPropagator."""
    uq = UncertaintyPropagator(sf=None)
    assert uq.sf is None
    assert isinstance(uq.model, UncertaintyModel)


@pytest.mark.fast
def test_uncertainty_propagator_add_conditions():
    """Test adding conditions to the propagator."""
    uq = UncertaintyPropagator(sf=None)
    uq.add_condition_uncertainty("Tgas", distribution="normal", std=50)
    uq.add_line_parameter_uncertainty(source="hitran")
    assert "Tgas" in uq.model.condition_params
    assert uq.model.has_line_uncertainty


@pytest.mark.fast
def test_uncertainty_propagator_propagate():
    """Test the Monte Carlo propagation loop with a mock SpectrumFactory."""
    from radis.test.tools.mock_sf import MockSpectrumFactory

    sf = MockSpectrumFactory()
    uq = UncertaintyPropagator(sf=sf)

    # Add a simple condition uncertainty
    uq.add_condition_uncertainty("Tgas", distribution="normal", std=10)

    # Run propagation with very few samples for speed
    result = uq.propagate(
        Tgas=1500,
        n_samples=5,
        confidence_level=0.95,
        method="monte_carlo",
        seed=42,
    )

    assert isinstance(result, ConfidenceBandResult)
    assert result.n_samples == 5
    assert result.quantity == "radiance_noslit"
    assert "mean" in result.get_stats()
    # verify w was extracted correctly
    assert len(result.wavenumber) == 100


@pytest.mark.fast
def test_uncertainty_propagator_parallel():
    """Test the Monte Carlo propagation loop with joblib parallelization."""
    pytest.importorskip("joblib", reason="joblib not installed")
    from radis.test.tools.mock_sf import MockSpectrumFactory

    sf = MockSpectrumFactory()
    uq = UncertaintyPropagator(sf=sf)

    # Add condition uncertainty (line uncertainty must be disabled for parallel)
    uq.add_condition_uncertainty("Tgas", distribution="normal", std=10)

    # Run propagation with n_jobs=2
    result = uq.propagate(
        Tgas=1500,
        n_samples=6,
        n_jobs=2,
        seed=42,
    )

    assert isinstance(result, ConfidenceBandResult)
    assert result.n_samples == 6


# ============================================================================
# Phase 2: Line-Parameter Perturbation Tests
# ============================================================================


@pytest.mark.fast
def test_line_parameter_perturbation():
    """Test the perturb/restore cycle on df0."""
    from radis.test.tools.mock_sf import MockSpectrumFactory

    sf = MockSpectrumFactory(n_lines=10)
    uq = UncertaintyPropagator(sf=sf)
    uq.add_line_parameter_uncertainty(source="hitran")

    # Save originals
    orig_int = sf.df0["int"].values.copy()
    orig_air = sf.df0["airbrd"].values.copy()

    rng = np.random.default_rng(123)
    uq._perturb_line_parameters(rng)

    # Values should be different after perturbation
    # Use array_equal (exact match) since perturbation adds noise
    assert not np.array_equal(
        sf.df0["int"].values, orig_int
    ), "int column should be perturbed"

    # Restore
    uq._restore_line_parameters()
    assert np.allclose(sf.df0["int"].values, orig_int), "int column should be restored"
    assert np.allclose(
        sf.df0["airbrd"].values, orig_air
    ), "airbrd column should be restored"


@pytest.mark.fast
def test_propagate_with_line_uncertainty():
    """Test MC propagation with line-parameter perturbation enabled."""
    from radis.test.tools.mock_sf import MockSpectrumFactory

    sf = MockSpectrumFactory()
    uq = UncertaintyPropagator(sf=sf)
    uq.add_line_parameter_uncertainty(source="hitran")
    uq.add_condition_uncertainty("Tgas", distribution="normal", std=10)

    result = uq.propagate(
        Tgas=1500,
        n_samples=5,
        seed=42,
    )

    assert isinstance(result, ConfidenceBandResult)
    assert result.n_samples == 5

    # After propagation df0 should be restored
    # (just verify no crash and result shape is correct)
    assert result.spectra_array.shape == (5, 100)


@pytest.mark.fast
def test_confidence_band_plot():
    """Test that plot_confidence_bands returns a matplotlib figure."""
    import matplotlib

    matplotlib.use("Agg")  # non-interactive backend

    rng = np.random.default_rng(42)
    wavenumber = np.linspace(2000, 2300, 50)
    base = np.sin(wavenumber / 100)
    spectra = base + rng.normal(0, 0.1, (20, 50))

    result = ConfidenceBandResult(
        wavenumber=wavenumber,
        spectra_array=spectra,
        quantity="radiance_noslit",
    )

    fig = result.plot_confidence_bands()
    assert fig is not None
    import matplotlib.pyplot as plt

    plt.close(fig)

    # Also test with show_individual
    fig2 = result.plot_confidence_bands(show_individual=True)
    assert fig2 is not None
    plt.close(fig2)


@pytest.mark.fast
def test_latin_hypercube_sampling():
    """Test that LHS produces correct sample shapes and coverage."""
    rng = np.random.default_rng(42)
    samples = UncertaintyPropagator._latin_hypercube_samples(
        n_samples=64,
        n_dims=3,
        rng=rng,
    )

    assert samples.shape == (64, 3)
    # All values should be in [0, 1)
    assert np.all(samples >= 0.0)
    assert np.all(samples < 1.0)

    # Each column should have good coverage of [0, 1]
    for d in range(3):
        col = samples[:, d]
        assert col.min() < 0.05, "LHS should cover near 0"
        assert col.max() > 0.95, "LHS should cover near 1"


@pytest.mark.fast
def test_propagate_latin_hypercube():
    """Test propagation with method='latin_hypercube'."""
    from radis.test.tools.mock_sf import MockSpectrumFactory

    sf = MockSpectrumFactory()
    uq = UncertaintyPropagator(sf=sf)

    # Test multiple distributions to hit all LHS code paths
    uq.add_condition_uncertainty("Tgas", distribution="normal", std=10)
    uq.add_condition_uncertainty(
        "mole_fraction", distribution="uniform", low=0.08, high=0.12
    )
    uq.add_condition_uncertainty("pressure", distribution="lognormal", s=0.1)

    result = uq.propagate(
        Tgas=1500,
        n_samples=8,
        method="latin_hypercube",
        seed=42,
    )

    assert isinstance(result, ConfidenceBandResult)
    assert result.n_samples == 8


# ============================================================================
# Phase 3: SensitivityAnalyzer Tests
# ============================================================================


@pytest.mark.fast
def test_sensitivity_analyzer_sobol_with_mock():
    """Test Sobol analysis with MockSpectrumFactory (requires SALib)."""
    pytest.importorskip("SALib", reason="SALib not installed")
    from radis.test.tools.mock_sf import MockSpectrumFactory

    sf = MockSpectrumFactory()
    sa = SensitivityAnalyzer(sf=sf)
    sa.add_parameter("Tgas", bounds=[1400, 1600])

    results = sa.sobol_analysis(n_samples=16)

    assert "S1" in results
    assert "ST" in results
    assert "wavenumber" in results
    assert results["S1"].shape[0] == 1  # 1 parameter
    assert results["parameter_names"] == ["Tgas"]


@pytest.mark.fast
def test_sensitivity_analyzer_error_budget():
    """Test error budget computation."""
    pytest.importorskip("SALib", reason="SALib not installed")
    from radis.test.tools.mock_sf import MockSpectrumFactory

    sf = MockSpectrumFactory()
    sa = SensitivityAnalyzer(sf=sf)
    sa.add_parameter("Tgas", bounds=[1400, 1600])

    results = sa.sobol_analysis(n_samples=16)
    budget = sa.get_error_budget(results)

    assert "Tgas" in budget
    # With a single parameter, its normalised contribution should be 1.0
    assert abs(budget["Tgas"] - 1.0) < 1e-10


@pytest.mark.fast
def test_sensitivity_analyzer_no_salib():
    """Test graceful error when SALib is not installed."""
    import sys
    import types

    # Temporarily block SALib
    saved = {}
    keys_to_block = [k for k in sys.modules if k == "SALib" or k.startswith("SALib.")]
    for k in keys_to_block:
        saved[k] = sys.modules.pop(k)

    blocker = types.ModuleType("SALib")
    blocker.__path__ = []

    def _raise(*a, **kw):
        raise ImportError("blocked")

    blocker.__getattr__ = _raise
    sys.modules["SALib"] = blocker
    sys.modules["SALib.analyze"] = None
    sys.modules["SALib.analyze.sobol"] = None
    sys.modules["SALib.sample"] = None
    sys.modules["SALib.sample.sobol"] = None

    try:
        sa = SensitivityAnalyzer(sf=None)
        sa.add_parameter("Tgas", bounds=[1200, 1800])
        with pytest.raises(ImportError, match="SALib"):
            sa.sobol_analysis()
    finally:
        for k in list(sys.modules.keys()):
            if k == "SALib" or k.startswith("SALib."):
                del sys.modules[k]
        sys.modules.update(saved)


@pytest.mark.fast
def test_sensitivity_analyzer_no_params():
    """Test that sobol_analysis with no parameters raises ValueError."""
    pytest.importorskip("SALib", reason="SALib not installed")

    sa = SensitivityAnalyzer(sf=None)
    with pytest.raises(ValueError, match="No parameters"):
        sa.sobol_analysis()


@pytest.mark.fast
def test_sensitivity_analyzer_plot():
    """Test that plot_sensitivity returns a matplotlib figure."""
    pytest.importorskip("SALib", reason="SALib not installed")
    import matplotlib

    matplotlib.use("Agg")

    import matplotlib.pyplot as plt

    from radis.test.tools.mock_sf import MockSpectrumFactory

    sf = MockSpectrumFactory()
    sa = SensitivityAnalyzer(sf=sf)
    sa.add_parameter("Tgas", bounds=[1400, 1600])

    results = sa.sobol_analysis(n_samples=16)

    # Test bar plot
    fig = sa.plot_sensitivity(results, kind="bar")
    assert fig is not None
    plt.close(fig)

    # Test heatmap plot
    fig_heat = sa.plot_sensitivity(results, kind="heatmap")
    assert fig_heat is not None
    plt.close(fig_heat)

    # Test unknown plot
    with pytest.raises(ValueError):
        sa.plot_sensitivity(results, kind="unknown")


if __name__ == "__main__":
    test_decode_hitran_uncertainty_all_codes()
    test_decode_hitran_uncertainty_position()
    test_decode_hitran_uncertainty_relative()
    test_decode_hitran_uncertainty_pressure_shift()
    test_decode_hitran_uncertainty_invalid_code()
    test_decode_hitran_uncertainty_invalid_parameter()
    test_decode_hitran_uncertainty_all_params()
    test_parse_hitran_ierr()
    test_parse_hitran_ierr_zero_padded()
    test_decode_geisa_uncertainty_valid()
    test_decode_geisa_uncertainty_zero()
    test_decode_geisa_uncertainty_nan()
    test_decode_geisa_uncertainty_none()
    test_uncertainty_model_creation()
    test_uncertainty_model_add_line_uncertainty()
    test_uncertainty_model_add_line_uncertainty_invalid()
    test_uncertainty_model_add_condition_normal()
    test_uncertainty_model_add_condition_uniform()
    test_uncertainty_model_add_condition_missing_params()
    test_uncertainty_model_sample_condition_normal()
    test_uncertainty_model_sample_condition_uniform()
    test_uncertainty_model_sample_condition_not_defined()
    test_uncertainty_model_invalid_distribution()
    test_confidence_band_result()
    test_spectrum_confidence_band_storage()
    test_spectrum_confidence_band_not_stored()
    test_sensitivity_analyzer_add_parameter()
    test_sensitivity_analyzer_invalid_bounds()
    test_get_auto_drop_columns_default()
    test_get_auto_drop_columns_preserve_uncertainty()
    test_get_auto_drop_columns_hitemp_preserve()
    test_get_auto_drop_columns_geisa()
    test_uncertainty_propagator_creation()
    test_uncertainty_propagator_add_conditions()
    test_uncertainty_propagator_propagate()
    test_line_parameter_perturbation()
    test_propagate_with_line_uncertainty()
    test_confidence_band_plot()
    test_latin_hypercube_sampling()
    test_propagate_latin_hypercube()
    print("All uncertainty tests passed!")
