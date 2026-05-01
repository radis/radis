# -*- coding: utf-8 -*-
"""
Spectral Uncertainty Quantification (UQ) Module
================================================

This module provides tools for propagating spectroscopic parameter
uncertainties through line-by-line simulations, producing spectra
with confidence bands and sensitivity analysis.

**Key components:**

- :func:`decode_hitran_uncertainty` — Decode HITRAN ``ierr`` codes
  into usable uncertainty ranges.
- :func:`decode_geisa_uncertainty` — Decode GEISA uncertainty columns
  into usable uncertainty ranges.
- :class:`UncertaintyModel` — Define parameter uncertainty distributions.
- :class:`UncertaintyPropagator` — Monte Carlo uncertainty propagation.
- :class:`SensitivityAnalyzer` — Sobol' global sensitivity analysis.
- :class:`ConfidenceBandResult` — Container for propagation results.

Notes
-----
No existing spectroscopy code provides built-in uncertainty propagation.
RADIS already parses uncertainty codes from HITRAN (``ierr``) and GEISA
(``ierrA``..``ierrU``), but drops them by default in
:data:`~radis.lbl.loader.drop_auto_columns_for_dbformat`. This module
turns that metadata into a scientifically impactful feature.

References
----------
.. [1] Gordon, I.E., et al. (2022). "The HITRAN2020 molecular
   spectroscopic database." JQSRT, 277, 107949.
.. [2] Pannier, E. & Laux, C.O. (2019). "RADIS: A nonequilibrium
   line-by-line radiative code." JQSRT, 222, 12-25.
.. [3] Sobol', I.M. (2001). "Global sensitivity indices for nonlinear
   mathematical models." Math. Comput. Simul., 55, 271-280.
"""

import numpy as np
from scipy import stats

# ============================================================================
# HITRAN Uncertainty Code Mapping
# ============================================================================

# For RELATIVE uncertainty parameters (intensity, air-broadened half-width,
# self-broadened half-width, temperature-dependence exponent):
#   code → (min_percent, max_percent)  where None means unbounded.
#
# For ABSOLUTE uncertainty parameters (line position, pressure shift):
#   code → (min_abs, max_abs) in cm⁻¹  where None means unbounded.
#
# Reference: Table 5 of Gordon et al. (2022), HITRAN2020.

HITRAN_IERR_RELATIVE = {
    0: None,  # Unreported or unavailable
    1: (20.0, None),  # Default or estimate, >= 20%
    2: (10.0, 20.0),  # Average or estimate
    3: (5.0, 10.0),
    4: (2.0, 5.0),
    5: (1.0, 2.0),
    6: (0.5, 1.0),
    7: (0.2, 0.5),
    8: (0.1, 0.2),
    9: (0.0, 0.1),  # Ideal
}
"""dict: HITRAN uncertainty codes for relative parameters
(intensity, half-widths, temperature-dependence).

Each code maps to ``(min_percent, max_percent)`` or ``None`` if
unreported. ``None`` in max_percent means unbounded (>= min).
"""

HITRAN_IERR_POSITION = {
    0: None,  # Unreported
    1: (1.0, None),  # >= 1 cm⁻¹
    2: (0.1, 1.0),  # 0.1 to 1.0 cm⁻¹
    3: (0.01, 0.1),
    4: (0.001, 0.01),
    5: (0.0001, 0.001),
    6: (0.00001, 0.0001),
    7: (0.000001, 0.00001),
    8: (0.0000001, 0.000001),
    9: (0.0, 0.0000001),  # Ideal
}
"""dict: HITRAN uncertainty codes for line position (absolute, cm⁻¹).

Each code maps to ``(min_cm1, max_cm1)`` or ``None`` if unreported.
"""

HITRAN_IERR_SHIFT = {
    0: None,  # Unreported
    1: (0.01, None),  # >= 0.01 cm⁻¹
    2: (0.001, 0.01),
    3: (0.0001, 0.001),
    4: (0.00001, 0.0001),
    5: (0.000001, 0.00001),
    6: (0.0000001, 0.000001),
    7: (0.0, 0.0000001),
    8: (0.0, 0.0000001),  # Same as 7 for shifts
    9: (0.0, 0.0000001),  # Ideal
}
"""dict: HITRAN uncertainty codes for pressure-induced line shift
(absolute, cm⁻¹.atm⁻¹)."""

# ierr digit index → parameter name mapping
HITRAN_IERR_PARAMS = {
    0: "position",  # ν (line position, cm⁻¹)
    1: "intensity",  # S (line intensity)
    2: "air_broadened",  # γ_air (air-broadened half-width)
    3: "self_broadened",  # γ_self (self-broadened half-width)
    4: "temperature_dependence",  # n_air (temperature exponent)
    5: "pressure_shift",  # δ_air (pressure-induced shift)
}
"""dict: Mapping from ierr digit index to parameter name."""


def decode_hitran_uncertainty(code, parameter="intensity"):
    """Decode a single HITRAN ``ierr`` uncertainty code digit.

    Parameters
    ----------
    code : int
        Uncertainty code (0–9). See Table 5 of [HITRAN2020]_.
    parameter : str
        Which parameter this code refers to. One of:

        - ``'position'`` — line position (absolute, cm⁻¹)
        - ``'intensity'`` — line intensity (relative, %)
        - ``'air_broadened'`` — air-broadened half-width (relative, %)
        - ``'self_broadened'`` — self-broadened half-width (relative, %)
        - ``'temperature_dependence'`` — temperature exponent (relative, %)
        - ``'pressure_shift'`` — pressure-induced shift (absolute, cm⁻¹)

    Returns
    -------
    tuple or None
        ``(min_uncertainty, max_uncertainty)`` or ``None`` if unreported.

        - For ``'position'``: absolute uncertainty in cm⁻¹
        - For ``'pressure_shift'``: absolute uncertainty in cm⁻¹.atm⁻¹
        - For all others: relative uncertainty in %

    Raises
    ------
    ValueError
        If ``code`` is not 0–9 or ``parameter`` is unknown.

    Examples
    --------
    ::

        >>> decode_hitran_uncertainty(5, 'intensity')
        (1.0, 2.0)
        >>> decode_hitran_uncertainty(4, 'position')
        (0.001, 0.01)
        >>> decode_hitran_uncertainty(0, 'intensity')  # unreported
        None

    References
    ----------
    .. [HITRAN2020] Gordon et al. (2022), JQSRT 277, 107949.
    """
    code = int(code)
    if code < 0 or code > 9:
        raise ValueError(f"HITRAN uncertainty code must be 0–9, got {code}")

    valid_params = set(HITRAN_IERR_PARAMS.values())
    if parameter not in valid_params:
        raise ValueError(
            f"Unknown parameter '{parameter}'. "
            f"Expected one of {sorted(valid_params)}"
        )

    if parameter == "position":
        return HITRAN_IERR_POSITION[code]
    elif parameter == "pressure_shift":
        return HITRAN_IERR_SHIFT[code]
    else:
        return HITRAN_IERR_RELATIVE[code]


def parse_hitran_ierr(ierr_value):
    """Parse a full HITRAN ``ierr`` 6-digit code into individual
    parameter uncertainties.

    Parameters
    ----------
    ierr_value : int or str
        The 6-digit ierr code. If int, zero-padded to 6 digits.

    Returns
    -------
    dict
        Dictionary mapping parameter names to ``(min, max)`` tuples.
        Parameters with unreported uncertainty (code 0) have value
        ``None``.

    Examples
    --------
    ::

        >>> parse_hitran_ierr(554400)
        {'position': (0.0001, 0.001),
         'intensity': (0.0001, 0.001),
         'air_broadened': (2.0, 5.0),
         'self_broadened': (2.0, 5.0),
         'temperature_dependence': None,
         'pressure_shift': None}
    """
    ierr_str = str(int(ierr_value)).zfill(6)
    result = {}
    for idx, param_name in HITRAN_IERR_PARAMS.items():
        digit = int(ierr_str[idx])
        result[param_name] = decode_hitran_uncertainty(digit, param_name)
    return result


# ============================================================================
# GEISA Uncertainty Decoder
# ============================================================================

# GEISA stores uncertainties as direct values (not codes)
GEISA_IERR_COLUMNS = {
    "ierrA": "position",  # line position uncertainty (cm⁻¹)
    "ierrB": "intensity",  # intensity uncertainty
    "ierrC": "air_broadened",  # air collision halfwidth unc.
    "ierrF": "temperature_dependence",  # temp-dep air halfwidth
    "ierrO": "pressure_shift",  # air pressure shift unc.
    "ierrR": "shift_temperature_dependence",  # temp-dep air shift
    "ierrN": "self_broadened",  # self-broadened unc.
    "ierrS": "self_temperature_dependence",  # temp-dep self halfwidth
    "ierrT": "self_pressure_shift",  # self pressure shift unc.
    "ierrU": "self_shift_temperature_dependence",  # temp-dep self shift
}
"""dict: Mapping from GEISA uncertainty column names to parameter names."""


def decode_geisa_uncertainty(ierr_value, parameter="position"):
    """Decode a GEISA uncertainty value.

    Unlike HITRAN (which uses integer codes 0–9), GEISA stores
    uncertainties as direct numerical values in dedicated columns.

    Parameters
    ----------
    ierr_value : float
        The uncertainty value from a GEISA ``ierrX`` column.
    parameter : str
        Which parameter this corresponds to (for context).

    Returns
    -------
    float or None
        The uncertainty value, or ``None`` if ``ierr_value`` is
        NaN or zero (indicating unreported).
    """
    if ierr_value is None or (isinstance(ierr_value, float) and np.isnan(ierr_value)):
        return None
    if ierr_value == 0.0:
        return None
    return float(ierr_value)


# ============================================================================
# UncertaintyModel
# ============================================================================


class UncertaintyModel:
    """Define uncertainty distributions for spectral parameters.

    Collects uncertainty information from database codes and/or
    user-specified distributions for experimental conditions.

    Parameters
    ----------
    source : str or None
        Default uncertainty source. ``'hitran'`` or ``'geisa'``.
        Default ``None``.

    Examples
    --------
    ::

        model = UncertaintyModel()
        model.add_line_uncertainty(source='hitran')
        model.add_condition('Tgas', distribution='normal', std=50)
        model.add_condition('pressure', distribution='uniform',
                            low=0.9, high=1.1)
    """

    def __init__(self, source=None):
        self.source = source
        self._line_uncertainty_source = None
        self._conditions = {}

    def add_line_uncertainty(self, source="hitran"):
        """Register that line parameter uncertainties should be
        read from database uncertainty codes.

        Parameters
        ----------
        source : str
            ``'hitran'`` or ``'geisa'``.
        """
        if source not in ("hitran", "geisa"):
            raise ValueError(
                f"Unknown uncertainty source '{source}'. "
                "Expected 'hitran' or 'geisa'."
            )
        self._line_uncertainty_source = source

    def add_condition(self, param, distribution="normal", **kwargs):
        """Add uncertainty on an experimental/simulation condition.

        Parameters
        ----------
        param : str
            Condition parameter name (e.g., ``'Tgas'``,
            ``'mole_fraction'``, ``'pressure'``, ``'path_length'``).
        distribution : str
            Distribution type. One of ``'normal'``, ``'uniform'``,
            ``'lognormal'``, ``'truncnorm'``.
        **kwargs
            Distribution parameters:

            - ``'normal'``: ``std`` (standard deviation)
            - ``'uniform'``: ``low``, ``high`` (bounds)
            - ``'lognormal'``: ``s`` (shape), ``scale``
            - ``'truncnorm'``: ``std``, ``low``, ``high``
        """
        dist = self._make_distribution(distribution, **kwargs)
        self._conditions[param] = {
            "distribution": distribution,
            "scipy_dist": dist,
            "kwargs": kwargs,
        }

    @property
    def condition_params(self):
        """list: Names of conditions with uncertainty."""
        return list(self._conditions.keys())

    @property
    def has_line_uncertainty(self):
        """bool: Whether line parameter uncertainty is enabled."""
        return self._line_uncertainty_source is not None

    @property
    def line_uncertainty_source(self):
        """str or None: Source of line parameter uncertainties."""
        return self._line_uncertainty_source

    def sample_condition(self, param, nominal_value, n_samples=1, seed=None):
        """Sample perturbed values for a condition parameter.

        Parameters
        ----------
        param : str
            Parameter name.
        nominal_value : float
            Nominal (central) value.
        n_samples : int
            Number of samples.
        seed : int or None
            Random seed.

        Returns
        -------
        numpy.ndarray
            Array of ``n_samples`` perturbed values.
        """
        if param not in self._conditions:
            raise KeyError(f"No uncertainty defined for condition '{param}'")

        rng = np.random.default_rng(seed)
        info = self._conditions[param]
        dist_type = info["distribution"]
        kwargs = info["kwargs"]

        if dist_type == "normal":
            std = kwargs["std"]
            return rng.normal(nominal_value, std, size=n_samples)
        elif dist_type == "uniform":
            low = kwargs["low"]
            high = kwargs["high"]
            return rng.uniform(low, high, size=n_samples)
        elif dist_type == "lognormal":
            s = kwargs["s"]
            scale = kwargs.get("scale", nominal_value)
            return rng.lognormal(np.log(scale), s, size=n_samples)
        elif dist_type == "truncnorm":
            std = kwargs["std"]
            low = kwargs.get("low", nominal_value - 4 * std)
            high = kwargs.get("high", nominal_value + 4 * std)
            a = (low - nominal_value) / std
            b = (high - nominal_value) / std
            return stats.truncnorm.rvs(
                a,
                b,
                loc=nominal_value,
                scale=std,
                size=n_samples,
                random_state=rng,
            )
        else:
            raise ValueError(f"Unknown distribution type '{dist_type}'")

    @staticmethod
    def _make_distribution(distribution, **kwargs):
        """Create a scipy.stats distribution object."""
        if distribution == "normal":
            if "std" not in kwargs:
                raise ValueError("Normal distribution requires 'std' parameter")
            return stats.norm(scale=kwargs["std"])
        elif distribution == "uniform":
            if "low" not in kwargs or "high" not in kwargs:
                raise ValueError(
                    "Uniform distribution requires 'low' and " "'high' parameters"
                )
            loc = kwargs["low"]
            scale = kwargs["high"] - kwargs["low"]
            return stats.uniform(loc=loc, scale=scale)
        elif distribution == "lognormal":
            if "s" not in kwargs:
                raise ValueError("Lognormal distribution requires 's' parameter")
            return stats.lognorm(s=kwargs["s"])
        elif distribution == "truncnorm":
            if "std" not in kwargs:
                raise ValueError("Truncated normal requires 'std' parameter")
            std = kwargs["std"]
            low = kwargs.get("low", -4 * std)
            high = kwargs.get("high", 4 * std)
            a = low / std
            b = high / std
            return stats.truncnorm(a, b, scale=std)
        else:
            raise ValueError(
                f"Unknown distribution '{distribution}'. "
                "Supported: 'normal', 'uniform', 'lognormal', "
                "'truncnorm'."
            )


# ============================================================================
# ConfidenceBandResult
# ============================================================================


class ConfidenceBandResult:
    """Container for uncertainty propagation results.

    Stores the ensemble of spectra from Monte Carlo sampling and
    provides methods to compute statistics and plot confidence bands.

    Parameters
    ----------
    wavenumber : numpy.ndarray
        Wavenumber grid (cm⁻¹).
    spectra_array : numpy.ndarray
        2-D array of shape ``(n_samples, len(wavenumber))``,
        containing the spectral quantity for each MC sample.
    quantity : str
        Name of the spectral quantity (e.g., ``'radiance_noslit'``).
    confidence_level : float
        Confidence level for bands (e.g., 0.95). Default 0.95.
    conditions : dict or None
        Conditions used for the calculation.
    """

    def __init__(
        self,
        wavenumber,
        spectra_array,
        quantity="radiance_noslit",
        confidence_level=0.95,
        conditions=None,
    ):
        self.wavenumber = np.asarray(wavenumber)
        self.spectra_array = np.asarray(spectra_array)
        self.quantity = quantity
        self.confidence_level = confidence_level
        self.conditions = conditions or {}
        self._n_samples = self.spectra_array.shape[0]

    @property
    def n_samples(self):
        """int: Number of Monte Carlo samples."""
        return self._n_samples

    @property
    def mean(self):
        """numpy.ndarray: Mean spectrum."""
        return np.mean(self.spectra_array, axis=0)

    @property
    def std(self):
        """numpy.ndarray: Standard deviation at each wavenumber."""
        return np.std(self.spectra_array, axis=0, ddof=1)

    @property
    def lower(self):
        """numpy.ndarray: Lower confidence band."""
        alpha = 1.0 - self.confidence_level
        return np.percentile(self.spectra_array, 100.0 * alpha / 2.0, axis=0)

    @property
    def upper(self):
        """numpy.ndarray: Upper confidence band."""
        alpha = 1.0 - self.confidence_level
        return np.percentile(self.spectra_array, 100.0 * (1.0 - alpha / 2.0), axis=0)

    def get_confidence_band(self):
        """Return the confidence band.

        Returns
        -------
        tuple
            ``(wavenumber, lower, upper)`` arrays.
        """
        return self.wavenumber, self.lower, self.upper

    def get_stats(self):
        """Return summary statistics at each wavenumber.

        Returns
        -------
        dict
            Keys: ``'mean'``, ``'std'``, ``'lower'``, ``'upper'``,
            ``'median'``, ``'n_samples'``, ``'confidence_level'``.
        """
        return {
            "mean": self.mean,
            "std": self.std,
            "lower": self.lower,
            "upper": self.upper,
            "median": np.median(self.spectra_array, axis=0),
            "n_samples": self.n_samples,
            "confidence_level": self.confidence_level,
        }


# ============================================================================
# UncertaintyPropagator
# ============================================================================


class UncertaintyPropagator:
    """Monte Carlo uncertainty propagation for spectral simulations.

    Wraps a :class:`~radis.lbl.factory.SpectrumFactory` and runs
    multiple simulations with perturbed parameters to estimate
    uncertainty in the output spectrum.

    Parameters
    ----------
    sf : SpectrumFactory
        A configured SpectrumFactory with a loaded databank.
    model : UncertaintyModel or None
        Pre-configured uncertainty model. If ``None``, a new one
        is created.

    Examples
    --------
    ::

        from radis import SpectrumFactory
        from radis.tools.uncertainty import (
            UncertaintyPropagator, UncertaintyModel
        )

        sf = SpectrumFactory(...)
        sf.fetch_databank('hitran', load_columns='all')

        uq = UncertaintyPropagator(sf)
        uq.model.add_condition('Tgas', distribution='normal', std=50)

        result = uq.propagate(Tgas=1500, n_samples=100)
        w, lower, upper = result.get_confidence_band()

    Notes
    -----
    Full Monte Carlo propagation (Phase 2) will iterate over
    ``n_samples``, perturbing parameters according to the
    :class:`UncertaintyModel` and collecting the resulting spectra.

    The current Phase 1 implementation provides the framework and
    supports condition-parameter perturbation. Line-parameter
    perturbation (perturbing individual ``S``, ``gamma_air``, etc.
    based on ``ierr`` codes) will be added in Phase 2.
    """

    def __init__(self, sf, model=None):
        self.sf = sf
        self.model = model or UncertaintyModel()

    def add_line_parameter_uncertainty(self, source="hitran"):
        """Enable line parameter uncertainty from database codes.

        Parameters
        ----------
        source : str
            ``'hitran'`` or ``'geisa'``.
        """
        self.model.add_line_uncertainty(source=source)

    def add_condition_uncertainty(self, param, distribution="normal", **kwargs):
        """Add uncertainty on a simulation condition.

        Parameters
        ----------
        param : str
            Condition name (e.g., ``'Tgas'``, ``'mole_fraction'``).
        distribution : str
            Distribution type.
        **kwargs
            Distribution parameters.
        """
        self.model.add_condition(param, distribution, **kwargs)

    def propagate(
        self,
        Tgas,
        n_samples=200,
        spectral_quantity="radiance_noslit",
        confidence_level=0.95,
        method="monte_carlo",
        seed=None,
        **eq_spectrum_kwargs,
    ):
        """Run uncertainty propagation via Monte Carlo sampling.

        Parameters
        ----------
        Tgas : float
            Nominal gas temperature (K).
        n_samples : int
            Number of Monte Carlo samples. Default 200.
        spectral_quantity : str
            Spectral quantity to propagate. Default
            ``'radiance_noslit'``.
        confidence_level : float
            Confidence level for bands. Default 0.95.
        method : str
            Sampling method: ``'monte_carlo'`` or
            ``'latin_hypercube'``. Default ``'monte_carlo'``.
        seed : int or None
            Random seed for reproducibility.
        **eq_spectrum_kwargs
            Additional keyword arguments passed to
            ``sf.eq_spectrum()``.

        Returns
        -------
        ConfidenceBandResult
            Results container with mean, std, and confidence bands.

        Notes
        -----
        Currently supports condition-parameter perturbation only.
        Line-parameter perturbation will be implemented in Phase 2.
        """
        rng = np.random.default_rng(seed)

        # Prepare nominal conditions
        nominal_conditions = {"Tgas": Tgas}
        nominal_conditions.update(eq_spectrum_kwargs)

        # Generate perturbed condition samples
        condition_samples = {}
        for param in self.model.condition_params:
            if param in nominal_conditions:
                nominal = nominal_conditions[param]
            elif param == "mole_fraction":
                nominal = self.sf.input.mole_fraction
            elif param == "pressure":
                nominal = self.sf.input.pressure
            elif param == "path_length":
                nominal = self.sf.input.path_length
            else:
                raise ValueError(
                    f"Cannot find nominal value for " f"parameter '{param}'"
                )
            condition_samples[param] = self.model.sample_condition(
                param,
                nominal,
                n_samples=n_samples,
                seed=rng.integers(0, 2**31),
            )

        # Run MC loop
        spectra_list = []
        wavenumber = None

        for i in range(n_samples):
            # Build kwargs for this sample
            sample_kwargs = dict(nominal_conditions)
            for param, samples in condition_samples.items():
                sample_kwargs[param] = float(samples[i])

            # Calculate spectrum
            s = self.sf.eq_spectrum(**sample_kwargs)
            w, I = s.get(spectral_quantity)

            if wavenumber is None:
                wavenumber = w
            spectra_list.append(I)

        spectra_array = np.array(spectra_list)

        return ConfidenceBandResult(
            wavenumber=wavenumber,
            spectra_array=spectra_array,
            quantity=spectral_quantity,
            confidence_level=confidence_level,
            conditions=nominal_conditions,
        )


# ============================================================================
# SensitivityAnalyzer (stub — full implementation in Phase 3)
# ============================================================================


class SensitivityAnalyzer:
    """Global sensitivity analysis for spectral simulations.

    Uses Sobol' indices (via SALib) to determine which input
    parameters contribute most to spectral uncertainty at each
    wavenumber.

    Parameters
    ----------
    sf : SpectrumFactory
        A configured SpectrumFactory with a loaded databank.

    Notes
    -----
    This is a Phase 3 stub. Full implementation requires the
    optional ``SALib`` dependency.

    Examples
    --------
    ::

        from radis.tools.uncertainty import SensitivityAnalyzer

        sa = SensitivityAnalyzer(sf)
        sa.add_parameter('Tgas', bounds=[1200, 1800])
        sa.add_parameter('mole_fraction', bounds=[0.05, 0.15])
        indices = sa.sobol_analysis(n_samples=1024)
    """

    def __init__(self, sf):
        self.sf = sf
        self._parameters = {}

    def add_parameter(self, name, bounds):
        """Add a parameter for sensitivity analysis.

        Parameters
        ----------
        name : str
            Parameter name (e.g., ``'Tgas'``).
        bounds : list
            ``[lower_bound, upper_bound]``.
        """
        if len(bounds) != 2:
            raise ValueError(f"bounds must be [low, high], got {bounds}")
        if bounds[0] >= bounds[1]:
            raise ValueError(f"Lower bound must be < upper bound, " f"got {bounds}")
        self._parameters[name] = {
            "bounds": list(bounds),
        }

    @property
    def parameter_names(self):
        """list: Names of parameters added for analysis."""
        return list(self._parameters.keys())

    def sobol_analysis(
        self,
        n_samples=1024,
        spectral_quantity="radiance_noslit",
    ):
        """Run Sobol' sensitivity analysis.

        Parameters
        ----------
        n_samples : int
            Number of Sobol' samples (must be power of 2).
        spectral_quantity : str
            Spectral quantity to analyze.

        Returns
        -------
        dict
            Placeholder. Full implementation in Phase 3.

        Raises
        ------
        NotImplementedError
            This is a Phase 3 stub.
        """
        raise NotImplementedError(
            "SensitivityAnalyzer.sobol_analysis() will be "
            "implemented in Phase 3. Requires the optional "
            "'SALib' dependency."
        )
