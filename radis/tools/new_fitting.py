# -*- coding: utf-8 -*-
"""

Module for New RADIS Fitting Interface

This module is to conduct benchmarking process on spectra (LTE or Non-LTE) with
temperatures (gas, vibrational, rotational, translational), mole fraction, cutoff,
as well as other characteristics as fitting parameters.

This module is implemented into RADIS from the Fitting Benchmarking repository, as
a 2022 GSOC project: https://github.com/TranHuuNhatHuy/RADIS-Spectrum-Fitting-Benchmark

This module aims to provide RADIS with a unified fitting interface which covers most of
current fitting cases encountered from fellow scientists and researchers, while being
simple and easy-to-use to RADIS newcomers, as you only need to input the parameters into
pre-defined Dictionary objects in a simplified script, just as filling a form. This also
serves as a bridge between RADIS codebase and a future fitting feature on RADIS App.

"""


import json
import time
from typing import Union

import numpy as np
from lmfit import Parameters, fit_report, minimize
from lmfit.minimizer import MinimizerResult

from radis import Spectrum, SpectrumFactory
from radis.spectrum import plot_diff
from radis.spectrum.compare import get_residual
from radis.tools.database import load_spec

# --------------------- MODELS / COST FUNCTIONS --------------------- #


def residual_LTE(params, conditions, s_data, sf, log, verbose):
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
        if True, will print out result of each fitting loop.

    Returns
    -------
    residual: float
        residuals of the two spectra, using RADIS's get_residual().

    Examples
    --------

    .. minigallery:: radis.tools.new_fitting.residual_LTE

    See Also
    --------
    :py:func:`~radis.tools.new_fitting.residual_NonLTE`
    :py:func:`~radis.tools.fitting.LTEModel`

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

    return residual


def residual_NonLTE(params, conditions, s_data, sf, log, verbose):
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
        if True, will print out result of each fitting loop.

    Returns
    -------
    residual: float
        residuals of the two spectra, using RADIS's get_residual().

    Examples
    --------

    .. minigallery:: radis.tools.new_fitting.residual_NonLTE

    See Also
    --------
    :py:func:`~radis.tools.new_fitting.residual_LTE`
    :py:func:`~radis.tools.fitting.radis.tools.fitting.Tvib12Tvib3Trot_NonLTEModel`

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

    return residual


# ----------------------------------------------------------------------- #


# --------------------- FITTING METHODS / FUNCTIONS --------------------- #


def get_conditions(
    input_conditions,
    verbose=True,
) -> Union[dict, Parameters]:
    """Get the ground-truth information by reading from input - either a JSON file
    or a dict parameter. Return a dict parameter containing fixed conditions for
    spectrum calculation, and an LMFIT.Parameters object containing fit parameters.

    Parameters
    ----------
    input_conditions : str or dict
        If str, the code assumes this is path to the JSON input file from the script
        calling it, for example:: "./test_dir/CO2_measured_spectrum_4-5um.json", and
        then parses the JSON data as a dict.
        If dict, the code just simply parses the dict, yeah, that simple.

    Other Parameters
    ----------------
    verbose: bool
        by default, True, print details about the fitting progress.

    Examples
    --------

    .. minigallery:: radis.tools.new_fitting.get_conditions

    """

    # LOG THE INPUT CONDITIONS (JSON OR DICT)

    if type(input_conditions) == str:  # Must be a JSON path

        with open(input_conditions, "r") as f:  # Load JSON as dict
            conditions = json.load(f)
            f.close()

    elif type(input_conditions) == dict:  # Must be just a dict

        conditions = input_conditions

    else:

        raise TypeError(
            f"Unrecognized input. Got type {type(input_conditions)}, expecting str or dict."
        )

    # READ "fit" KEY AND LOAD IT INTO NEWLY INITIALIZED PARAMETERS OBJECT

    params = Parameters()  # Initiate Parameters object

    fit_params = conditions.pop("fit")
    dhb = 1000  # Default half-bound for fit temperatures, in case user don't specify bounds

    # Deal with "Tvib", especially if it has multiple temperatures

    if "Tvib" in fit_params:

        Tvib = fit_params["Tvib"]

        if (type(Tvib) == list) or (
            type(Tvib) == tuple
        ):  # Multiple temperatures detected

            # Split Tvib[] into Tvib0, Tvib1, Tvib2,... for later being put into Parameters
            for i in range(len(Tvib)):
                fit_params[f"Tvib{i}"] = Tvib[i]

            # Destroy the original "Tvib" key
            fit_params.pop("Tvib")

    # Deal with "offset", which is something like "-0.2 nm", as requested by Mr. Minou

    if "offset" in fit_params:

        offset_value, offset_unit = fit_params.pop("offset").split()
        if offset_unit == "cm-1":
            offset_unitName = "cm1"
        elif offset_unit == "nm":
            offset_unitName = "nm"
        else:
            raise ValueError(
                f'Unrecognized offset unit {offset_unit}. Only "nm" and "cm-1" are accepted.'
            )

        # The reason for all the mess above is, LMFIT does not accept parameter name with hyphen like cm-1.
        # Thus, I remove the hyphen when assigning its name into LMFIT Parameters object.

        # Get numerical value only, and add unit to name for later identification
        fit_params[f"offset{offset_unitName}"] = float(offset_value)

    # Put every fit parameter into Parameters, with the default bounding range

    for param in fit_params:

        init_val = fit_params[param]

        if param == "mole_fraction":
            init_bound = [0, 1]
        elif "offset" in param:  # Either "offsetnm" or "offsetcm1"
            init_bound = [-1, 1]
        else:
            init_bound = [(init_val - dhb) if init_val > dhb else 0, init_val + dhb]

        params.add(
            param,
            value=init_val,
            min=init_bound[0],
            max=init_bound[1],
        )

    # READ "bounds" KEY AND RE-ASSIGN DEFAULT BOUNDING RANGES WITH USER-DEFINED ONES

    if "bounds" in conditions:

        bounds = conditions.pop("bounds")

        # Deal with "Tvib", especially if it has multiple temperature-bounding ranges

        if "Tvib" in bounds:

            bound_Tvib = bounds["Tvib"]

            if (
                len(bound_Tvib) > 1
            ):  # Not sure if it's multiple ranges, or just 1 range with 2 ends
                if type(bound_Tvib[0]) in [
                    list,
                    tuple,
                ]:  # Now we sure it's multiple ranges

                    # Split Tvib[] into Tvib0, Tvib1, Tvib2,... for later being put into Parameters
                    for i in range(len(bound_Tvib)):
                        bounds[f"Tvib{i}"] = bound_Tvib[i]

                    # Destroy the original key
                    bounds.pop("Tvib")

        # Deal with "offset", especially now the offset parameter is either "offsetnm" or "offsetcm1"

        if "offset" in bounds:

            bound_offset = bounds.pop("offset")

            if "offsetnm" in params:
                bounds["offsetnm"] = bound_offset

            if "offsetcm1" in params:
                bounds["offsetcm1"] = bound_offset

        # Set user-defined bounding ranges

        for bound_param in bounds:

            bound_val = bounds[bound_param]

            params[bound_param].set(min=bound_val[0], max=bound_val[1])

    return conditions, params


# -------------------------------------------------------------------------------------------- #


def spectrum_refinement(s_data, conditions, verbose) -> Union[Spectrum, dict]:
    """Receive an experimental spectrum and further refine it according to the provided pipeline
    and ground-truths. Refinement process includes extracting the desired spectrum quantity,
    removing NaN values, and other additional convolutions. Finally, a refined Spectrum object is
    returned, along with the original condition as a Dictionary object.

    Parameters
    ----------
    s_data : Spectrum
        experimental spectrum loaded from the path given in input JSON file.
    conditions : dict
        a Dictionary containing all ground-truth information, and also spectrum refinements or
        convolutions (desired spectrum quantity, normalization, etc.), as well as fitting process
        (max loop allowed, terminal tolerance, etc.).

    Other parameters
    ----------------
    verbose : bool
        by default, True, print details about the refinement progress.

    Returns
    -------
    s_refined : Spectrum
        the refined spectrum.
    conditions : dict
        the input conditions that might have been added fit_var (in case user didn't).

    Examples
    --------

    .. minigallery:: radis.tools.new_fitting.spectrum_refinement

    """

    # Extract spectrum

    pipeline = conditions["pipeline"]

    if "fit_var" in pipeline:
        fit_var = pipeline["fit_var"]  # Acquire the stated quantity
    else:
        fit_var = s_data.get_vars()[0]  # If not stated, take first quantity
        conditions["pipeline"]["fit_var"] = fit_var  # And add to the dict

    s_data_wav, s_data_val = s_data.get(fit_var)

    if verbose:
        print(f"\nAcquired spectral quantity '{fit_var}' from the spectrum.")

    # Remove NaN values. A wise man once said, "Nan is good but only in India"

    s_data_mtr = np.vstack(
        (s_data_wav, s_data_val)
    )  # Merge wavelength/wavenumber and intensity into pairs
    s_data_mtr = s_data_mtr[:, ~np.isnan(s_data_mtr).any(axis=0)]  # Purge NaN pairs

    if verbose:
        print(f"NaN values successfully purged.", end="")
        print(f"Number of data points left: {len(s_data_mtr[0])} points.")

    # Recreate the data spectrum with the spectral quantity

    s_refined = (
        Spectrum(
            {fit_var: (s_data_mtr[0], s_data_mtr[1])},
            wunit=s_data.get_waveunit(),
            units={fit_var: s_data.units[fit_var]},
        )
        .take(fit_var)
        .sort()
    )

    # Further refinement

    # Apply normalization if required
    if "normalize" in pipeline:
        if pipeline["normalize"]:
            s_refined = s_refined.normalize()
            if verbose:
                print("Normalization applied.")

    # Currently supporting normalization. May add more in the future if requested.

    return s_refined, conditions


# -------------------------------------------------------------------------------------------- #


def fit_spectrum(
    s_exp=None,
    fit_params=None,
    model=None,
    pipeline=None,
    bounds=None,
    input_file=None,
    verbose=True,
    show_plot=True,
    fit_kws={},
) -> Union[Spectrum, MinimizerResult, dict]:
    """Fit an experimental spectrum (from here referred as "data spectrum") with a modeled one,
    then derive the fit results. Data spectrum is loaded from the path stated in the JSON file,
    while model spectrum is generated based on the conditions stated in the JSON file, too.

    Parameters
    ----------
    s_exp : Spectrum
        a Spectrum object, containing the experimental spectrum.
    fit_params : dict
        a dict object, containing fit parameters and their initial values.
    bounds : dict
        a dict object, optional, containing bounding ranges of those fit parameters. If not stated,
        all fit parameters will receive default bounding ranges.
    model : dict
        a dict object, containing experimental conditions and some additional convolutions for
        modeling the spectrum.
    pipeline : dict
        a dict object, containing some preference properties of fitting process.
    input_file : str
        path to the JSON input file from the script calling it, by default, None. If this parameter
        is defined, it means user uses JSON file for inputing everything, and so will discard all
        other input parameters. For example:: "./test_dir/CO2_measured_spectrum_4-5um.json"

    Other parameters
    ----------
    verbose : bool
        print details about each loop of the fitting process. By default, False.
    show_plot : bool
        show experimental spectrum before and after being fitted. By default, True. Showing the plot
        will pause the runtime, so if you are running benchmark which involves multiple fitting cases,
        set it to False.
    fit_kws: dict
        Options to pass to the minimizer being used. See :py:func:`~lmfit.minimizer.minimize`

    Returns
    -------
    s_best: Spectrum
        visualization of the best fit results obtained, as a spectrum.
    best_fit: MinimizerResult
        best fit results, output of LMFIT MinimizeResult.
    log: dict
        a Dictionary storing runtime log of the fitting process that are
        not quite covered by the Minimizer, including: residual and fit
        values after each fitting loop, and total time elapsed.

    Examples
    --------

    ::
        from radis import load_spec
        from radis.test.utils import getTestFile
        from radis.tools.new_fitting import fit_spectrum

        # Load an experimental spectrum. You can prepare yours, or fetch one of them in the radis/test/files directory.
        my_spec = getTestFile("synth-NH3-1-500-2000cm-P10-mf0.01-p1.spec")
        s_experimental = load_spec(my_spec)

        # Experimental conditions which will be used for spectrum modeling. Basically, these are known ground-truths.
        experimental_conditions = {
            "molecule": "NH3",  # Molecule ID
            "isotope": "1",  # Isotope ID, can have multiple at once
            "wmin": 1000,  # Starting wavelength/wavenumber to be cropped out from the original experimental spectrum.
            "wmax": 1050,  # Ending wavelength/wavenumber for the cropping range.
            "wunit": "cm-1",  # Accompanying unit of those 2 wavelengths/wavenumbers above.
            "mole_fraction": 0.01,  # Species mole fraction, from 0 to 1.
            "pressure": 10,  # Total pressure of gas, in "bar" unit by default, but you can also use Astropy units.
            "path_length": 1,  # Experimental path length, in "cm" unit by default, but you can also use Astropy units.
            "slit": "1 nm",  # Experimental slit, must be a blank space separating slit amount and unit.
            "offset": "-0.2 nm",  # Experimental offset, must be a blank space separating offset amount and unit.
            "wstep": 0.003,  # Resolution of wavenumber grid, in cm-1.
            "databank": "hitran",  # Databank used for the spectrum calculation. Must be stated.
        }

        # List of parameters to be fitted, accompanied by their initial values.
        fit_parameters = {
            "Tgas": 700,  # Gas temperature, in K.
        }

        # List of bounding ranges applied for those fit parameters above.
        bounding_ranges = {
            "Tgas": [
                500,
                2000,
            ],
        }

        # Fitting pipeline setups.
        fit_properties = {
            "method": "leastsq",  # Preferred fitting method. By default, "leastsq".
            "fit_var": "radiance",  # Spectral quantity to be extracted for fitting process, such as "radiance", "absorbance", etc.
            "normalize": False,  # Either applying normalization on both spectra or not.
            "max_loop": 150,  # Max number of loops allowed. By default, 200.
            "tol": 1e-15,  # Fitting tolerance, only applicable for "lbfgsb" method.
        }

        # Conduct the fitting process!
        s_best, result, log = fit_spectrum(
            s_exp=s_experimental,  # Experimental spectrum.
            fit_params=fit_parameters,  # Fit parameters.
            bounds=bounding_ranges,  # Bounding ranges for those fit parameters.
            model=experimental_conditions,  # Experimental ground-truths conditions.
            pipeline=fit_properties,  # Fitting pipeline references.
            verbose=False,  # If you want a clean result, stay False. If you want to see more about each loop, True.
        )
    .. minigallery:: radis.tools.new_fitting.fit_spectrum
    Notes
    -----
    For the fitting method, you can try one among 17 different fitting methods and algorithms of LMFIT, introduced
    in `LMFIT method list <https://lmfit.github.io/lmfit-py/fitting.html#choosing-different-fitting-methods>`.
    See Also
    --------
    :py:func:`~radis.tools.fitting.fit_spectrum()`,
    :py:meth:`~radis.lbl.factory.SpectrumFactory.fit_spectrum`,
    :py:mod:`fitroom`

    """

    begin = time.time()  # Start the fitting time counter

    if verbose:
        print(
            "\n======================= COMMENCE FITTING PROCESS ======================="
        )

    # ACQUIRE AND REFINE EXPERIMENTAL SPECTRUM s_data

    if input_file != None:  # The user uses JSON file

        # Load data from JSON file, then create a Parameters object
        conditions, params = get_conditions(input_file)

        # Get s_data spectrum from the path stated in acquired JSON data, assuming the JSON
        # and the .spec files are in the same directory (it SHOULD be!)

        fileName = conditions["model"].pop("fileName")
        if fileName[-5:] != ".spec":  # fileName doesn't have ".spec" extension
            raise ValueError('Warning: fileName must include ".spec" extension!')

        input_dir = "/".join(input_file.split("/")[0:-1]) + "/"
        spec_path = input_dir + fileName

        # Load and crop the experimental spectrum
        s_data = load_spec(spec_path).crop(
            conditions["model"]["wmin"],
            conditions["model"]["wmax"],
            conditions["model"]["wunit"],
        )

    else:  # The user states a bunch of dicts separately

        # Merge all separated dict inputs into a single dict first

        conditions = {
            "model": model,
            "fit": fit_params,
            "bounds": bounds if bounds is not None else {},
            "pipeline": pipeline,
        }

        # Load data, create Parameters object, from the above dict
        conditions, params = get_conditions(conditions)

        # Load and crop the experimental spectrum
        if "wunit" in conditions["model"]:  # Users state wunit directly
            s_data = s_exp.crop(
                conditions["model"]["wmin"],
                conditions["model"]["wmax"],
                conditions["model"]["wunit"],
            )
        else:  # Users don't state wunit directly, probably astropy units
            s_data = s_exp.crop(
                conditions["model"]["wmin"],
                conditions["model"]["wmax"],
            )

    # Log the time mark when experimental data is retrieved
    end_exp_load = time.time()
    time_exp_load = end_exp_load - begin

    if verbose:
        print(f"\nSuccessfully retrieved the experimental data in {time_exp_load}s.")

    # Further refine the data spectrum before calculating diff

    s_data, conditions = spectrum_refinement(s_data, conditions, verbose)

    pipeline = conditions["pipeline"]
    fit_var = pipeline["fit_var"]

    # Log the time mark when spectrum refinement is over
    end_exp_refine = time.time()
    time_exp_refine = end_exp_refine - end_exp_load

    if verbose:
        print(f"Successfully refined the experimental data in {time_exp_refine}s.")
    if show_plot:
        s_data.plot(show=True)

    # PRE-MINIMIZATION SETUP

    # Create a Spectrum Factory object for modeling spectra

    kwargs = {}
    model = conditions["model"]

    # List of keys that are included in "model", but not included in SpectrumFactory syntax
    ignore_keys = [
        "slit",
        "offset",
        "databank",
    ]

    for cond in model:
        if cond not in ignore_keys:
            kwargs[cond] = model[cond]

    sf = SpectrumFactory(
        **kwargs,
        verbose=False,
        warnings={
            "MissingSelfBroadeningWarning": "ignore",
            "UserWarning": "ignore",
            "NegativeEnergiesWarning": "ignore",
            "HighTemperatureWarning": "ignore",
        },
    )

    # Decide the type of model - LTE or Non-LTE
    LTE = True  # LTE == True means it's LTE
    if ("Tvib" in params) or ("Tvib0" in params):
        LTE = False  # LTE == False means it's non-LTE

    # Determine fitting method, either stated by user or "lbfgsb" by default
    if "method" in pipeline:
        method = pipeline["method"]  # User-defined
    else:
        method = "leastsq"  # By default

    # Determine additional fitting conditions for the minimizer

    fit_kws.update(
        {
            "max_nfev": pipeline["max_loop"]
            if "max_loop" in pipeline
            else 200,  # Default max number of loops: 200
        }
    )

    if "tol" in pipeline:
        if pipeline["method"] == "lbfgsb":
            fit_kws["tol"] = pipeline["tol"]
        else:
            print('"tol" parameter spotted but "method" is not "lbfgsb"!')

    # Prepare fitting log
    log = {"residual": [], "fit_vals": [], "time_fitting": 0}

    # COMMENCE THE FITTING PROCESS

    databank = model.pop("databank")

    # Log the time mark when fitting begins
    begin_fitting = time.time()

    # For LTE spectra
    if LTE:

        sf.fetch_databank(databank, load_columns="equilibrium")

        if verbose:
            print("\nCommence fitting process for LTE spectrum!\n")

        result = minimize(
            residual_LTE,
            params,
            method=method,
            args=[conditions, s_data, sf, log, verbose],
            **fit_kws,
        )

    # For non-LTE spectra
    if not (LTE):

        sf.fetch_databank(databank, load_columns="noneq")

        if verbose:
            print("\nCommence fitting process for non-LTE spectrum!")

        result = minimize(
            residual_NonLTE,
            params,
            method=method,
            args=[conditions, s_data, sf, log, verbose],
            **fit_kws,
        )

    # POST-MINIMIZATION REPORT

    print(fit_report(result), end="\n")  # Report the fitting result

    # REGENERATE THE BEST-FIT SPECTRUM, AS A RESULT FROM THE MINIMIZER

    # Load initial values of fit parameters
    fit_show = {}
    offset_unit = ""
    for name, param in result.params.items():
        if "offset" in name:  # "offset" detected
            offset_value = float(param.value)
            offset_unit = name[6:]
            # Deal with "cm1"
            if offset_unit == "cm1":
                offset_unit = "cm-1"
        else:  # Other parameters not related to offset
            fit_show[name] = float(param.value)
    fit_show["name"] = "best_fit"

    # Generate best fitted spectrum result
    if LTE:
        s_result = sf.eq_spectrum(**fit_show)
    else:
        s_result = sf.non_eq_spectrum(**fit_show)

    # In case offset is detected in parameters
    if offset_unit != "":
        s_result = s_result.offset(offset_value, offset_unit)

    # Apply slit stated in "model"
    if "slit" in model:

        slit_model = model["slit"]

        # The user uses simple format of "[value] [unit]", such as "-0.2 nm"
        if isinstance(slit_model, str):
            slit_val, slit_unit = slit_model.split()
            s_result = s_result.apply_slit(float(slit_val), slit_unit)

        # The user uses a dict with complex format of slit function, unit, shape, center wavespace, dispersion, etc.
        if isinstance(slit_model, dict):
            kwargs = {}
            for cond in slit_model:
                kwargs[cond] = slit_model[cond]
            s_result = s_result.apply_slit(**kwargs)

    # Take spectral quantity
    s_result = s_result.take(fit_var)

    # Apply offset stated in "model"
    if "offset" in model:
        off_val, off_unit = model["offset"].split()
        s_result = s_result.offset(float(off_val), off_unit)

    # Apply normalization
    if "normalize" in pipeline:
        if pipeline["normalize"]:
            s_result = s_result.normalize()

    # Log the time mark when fitting is over
    end_fitting = time.time()
    time_fitting = float(end_fitting - begin_fitting)

    if verbose:
        print(f"\nSuccessfully finished the fitting process in {time_fitting}s.")

    log["time_fitting"] = time_fitting

    # PLOT THE DIFFERENCE BETWEEN THE TWO

    # If show_plot = True by default, show the diff plot
    if show_plot:
        plot_diff(
            s_data,
            s_result,
            fit_var,
            method=["diff", "ratio"],
            show=True,
        )

    if verbose:
        print(
            "\n======================== END OF FITTING PROCESS ========================\n"
        )

    return s_result, result, log
