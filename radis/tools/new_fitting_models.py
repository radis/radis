# -*- coding: utf-8 -*-

from radis import get_residual

# General model for LTE spectra. This model will be heavily modified during
# benchmarking process to see which one works best.


def residual_LTE(params, conditions, s_data, sf, log, verbose=True):
    """A cost function that calculates an LTE spectrum based on the
    initial conditions and values of fit parameters, then returning a
    scalar containing difference between experimental and data spectra.

    Parameters
    ----------
    params : LMFIT.Parameters
        basically a dict containing fit parameters that will be varied
        during the minimization.
    conditions: dict
        a dict containing conditions (fixed parameters) for generating
        model spectra.
    s_data: Spectrum
        data spectrum, loaded from the directory stated in JSON file.
    sf: SpectrumFactory
        the SpectrumFactory object used for modeling spectra, generated
        before the minimize loop occurs.
    log: list
        a Dictionary storing runtime log of the fitting process that are
        not quite covered by the Minimizer, including: residual and fit
        values after each fitting loop, and total time elapsed.

    Other parameters
    ----------------
    verbose : bool
        by default, True, will print out information of fitting process.

    Returns
    -------
    residual: float
        residuals of the two spectra, using RADIS's get_residual().

    """

    # GENERATE LTE SPECTRUM BASED ON THOSE PARAMETERS

    # Load initial values of fit parameters

    # These keys must be ignored when parsing parameters. They will be dealt with separately
    ignore_keys = [
        "offsetnm",
        "offsetcm1",
    ]

    # Retrieve current values of parameters, for calculating model spectrum
    kwargs = {}
    for param in params:
        if param not in ignore_keys:
            kwargs[param] = float(params[param])

    # Model spectrum calculation
    s_model = sf.eq_spectrum(**kwargs)

    # Deal with "offset"
    if "offsetnm" in params:
        offset_value = float(params["offsetnm"])
        s_model = s_model.offset(offset_value, "nm")
    if "offsetcm1" in params:
        offset_value = float(params["offsetcm1"])
        s_model = s_model.offset(offset_value, "cm-1")

    # FURTHER REFINE THE MODELED SPECTRUM BEFORE CALCULATING DIFF

    pipeline = conditions["pipeline"]
    model = conditions["model"]

    # Apply slit stated in "model"
    if "slit" in model:

        slit_model = model["slit"]

        # The user uses simple format of "[value] [unit]", such as "-0.2 nm"
        if isinstance(slit_model, str):
            slit_val, slit_unit = slit_model.split()
            s_model = s_model.apply_slit(float(slit_val), slit_unit)

        # The user uses a dict with complex format of slit function, unit, shape, center wavespace, dispersion, etc.
        if isinstance(slit_model, dict):
            kwargs = {}
            for cond in slit_model:
                kwargs[cond] = slit_model[cond]
            s_model = s_model.apply_slit(**kwargs)

    # Take spectral quantity
    fit_var = pipeline["fit_var"]
    s_model = s_model.take(fit_var)

    # Apply offset stated in "model"
    if "offset" in model:
        off_val, off_unit = model["offset"].split()
        s_model = s_model.offset(float(off_val), off_unit)

    # Apply normalization
    if "normalize" in pipeline:
        if pipeline["normalize"]:
            s_model = s_model.normalize()

    # ACQUIRE AND RETURN DIFF, ALSO LOG FITTING HISTORY

    # Acquire diff
    residual = get_residual(s_model, s_data, fit_var, norm="L2", ignore_nan="True")

    # Log the current residual
    log["residual"].append(residual)

    # Log the current fitting values of fit parameters
    current_fitvals = []
    for param in params:
        current_fitvals.append(float(params[param]))
    log["fit_vals"].append(current_fitvals)

    # Print information of fitting process
    if verbose:
        for param in params:
            print(f"{param} = {float(params[param])}")
        print(f"\nResidual = {residual}\n")

    return residual


def residual_NonLTE(params, conditions, s_data, sf, log, verbose=True):
    """A cost function that calculates a non-LTE spectrum based on the
    initial conditions and values of fit parameters, then returning a
    scalar containing difference between experimental and data spectra.

    Parameters
    ----------
    params : LMFIT.Parameters
        basically a dict containing fit parameters that will be varied
        during the minimization.
    conditions: dict
        a dict containing conditions (fixed parameters) for generating
        model spectra.
    s_data: Spectrum
        data spectrum, loaded from the directory stated in JSON file.
    sf: SpectrumFactory
        the SpectrumFactory object used for modeling spectra, generated
        before the minimize loop occurs.
    log: list
        a Dictionary storing runtime log of the fitting process that are
        not quite covered by the Minimizer, including: residual and fit
        values after each fitting loop, and total time elapsed.

    Other parameters
    ----------------
    verbose : bool
        by default, True, will print out information of fitting process.

    Returns
    -------
    residual: float
        residuals of the two spectra, using RADIS's get_residual().

    """

    # GENERATE NON-LTE SPECTRUM BASED ON THOSE PARAMETERS

    # Load initial values of fit parameters

    # These keys must be ignored when parsing parameters. They will be dealt with separately
    ignore_keys = [
        "offsetnm",
        "offsetcm-1",
    ]

    # Retrieve current values of parameters, for calculating model spectrum
    kwargs = {}
    for param in params:
        if param not in ignore_keys:
            kwargs[param] = float(params[param])

    # Deal with the case of multiple Tvib temperatures
    if "Tvib0" in kwargs:  # There is trace of a "Tvib fragmentation" before
        Tvib = []
        for kw in kwargs:
            if "Tvib" in kw:  # Such as "Tvib0", "Tvib1" or so
                Tvib.append(kwargs[kw])  # Bring them altogether, uwu
                kwargs.pop(kw)  # Dispose the fragmented one in kwargs
        kwargs["Tvib"] = tuple(Tvib)  # Finally, we have the tuple of Tvib

    # Model spectrum calculation
    s_model = sf.non_eq_spectrum(**kwargs)

    # Deal with "offset"
    if "offsetnm" in params:
        offset_value = float(params["offsetnm"])
        s_model = s_model.offset(offset_value, "nm")
    if "offsetcm-1" in params:
        offset_value = float(params["offsetcm-1"])
        s_model = s_model.offset(offset_value, "cm-1")

    # FURTHER REFINE THE MODELED SPECTRUM BEFORE CALCULATING DIFF

    pipeline = conditions["pipeline"]
    model = conditions["model"]

    # Apply slit stated in "model"
    if "slit" in model:

        slit_model = model["slit"]

        # The user uses simple format of "[value] [unit]", such as "-0.2 nm"
        if isinstance(slit_model, str):
            slit_val, slit_unit = slit_model.split()
            s_model = s_model.apply_slit(float(slit_val), slit_unit)

        # The user uses a dict with complex format of slit function, unit, shape, center wavespace, dispersion, etc.
        if isinstance(slit_model, dict):
            kwargs = {}
            for cond in slit_model:
                kwargs[cond] = slit_model[cond]
            s_model = s_model.apply_slit(**kwargs)

    # Take spectral quantity
    fit_var = pipeline["fit_var"]
    s_model = s_model.take(fit_var)

    # Apply offset
    if "offset" in model:
        off_val, off_unit = model["offset"].split()
        s_model = s_model.offset(float(off_val), off_unit)

    # Apply normalization
    if "normalize" in pipeline:
        if pipeline["normalize"]:
            s_model = s_model.normalize()

    # ACQUIRE AND RETURN DIFF, ALSO LOG FITTING HISTORY

    # Acquire diff
    residual = get_residual(s_data, s_model, fit_var, norm="L2", ignore_nan="True")

    # Log the current residual
    log["residual"].append(residual)

    # Log the current fitting values of fit parameters
    current_fitvals = []
    for param in params:
        current_fitvals.append(float(params[param]))
    log["fit_vals"].append(current_fitvals)

    # Print information of fitting process
    if verbose:
        for param in params:
            print(f"{param} = {float(params[param])}")
        print(f"\nResidual = {residual}\n")

    return residual
