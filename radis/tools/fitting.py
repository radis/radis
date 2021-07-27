# -*- coding: utf-8 -*-
"""

Fitting Tools embedded RADIS, for simple cases

For more advanced cases, use Fitroom : https://github.com/radis/fitroom

https://user-images.githubusercontent.com/16088743/120166810-4ac14980-c1fd-11eb-9dd5-8fb037db8793.mp4


Originally in radis-examples : https://github.com/radis/radis-examples

"""


from os.path import join

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

from radis import Spectrum, SpectrumFactory
from radis.spectrum import plot_diff
from radis.spectrum.compare import get_diff, get_residual
from radis.test.utils import getValidationCase, setup_test_line_databases


# Calculate a new spectrum for given parameters:
def Tvib12TrotModel(factory, model_input):
    """A model returning a single-slab non-LTE spectrum with Tvib=(T12, T12, T3), Trot

    Parameters
    ----------
    model_input: dict
        input dictionary (typically: temperatures). Example::
            {'T12':value,
             'T3':value,
             'Trot':value}

    Returns
    -------

    s: :class:`~radis.spectrum.spectrum.Spectrum`
        calculated Spectrum

    Examples
    --------

    .. minigallery:: radis.tools.fitting.Tvib12TrotModel

    """

    # ... input should remain a dict
    T12 = model_input["T12"]
    T3 = model_input["T3"]
    Trot = model_input["Trot"]

    # ... create whatever model below (can have several slabs with SerialSlabs
    # ... or MergeSlabs, etc.)

    # >>> This is where the RADIS calculation is done!

    s = factory.non_eq_spectrum(
        (T12, T12, T3),
        Trot,
        Ttrans=Trot,
        vib_distribution="treanor",
        name="treanor. fit",
    )

    # <<<

    # ... output should be a Spectrum object
    return s


ite = 2


def fit_spectrum(
    factory, s_exp, model, fit_parameters, bounds={}, plot=True, maxiter=300
):
    """Fit an experimental spectrum with an arbitrary model and an arbitrary
    number of fit parameters.

    Parameters
    ----------
    s_exp : Spectrum
        experimental spectrum. Should have only spectral array only. Use
        :py:meth:`~radis.spectrum.spectrum.Spectrum.take`, e.g::
            sf.fit_spectrum(s_exp.take('transmittance'))
    model : func -> Spectrum
        a line-of-sight model returning a Spectrum. Example : :py:func:`~radis.tools.fitting.Tvib12TrotModel`
    fit_parameters : dict
        ::
            {fit_parameter:initial_value}
    bounds : dict, optional
        ::
            {fit_parameter:[min, max]}

    Other Parameters
    ----------------
    plot: bool
        if True, plot spectra as they are computed; and plot the convergence of
        the residual.
    maxiter: int
        max number of iteration, default 300

    Returns
    -------
    s_best

    See Also
    --------
    :py:meth:`~radis.lbl.factory.SpectrumFactory.fit_spectrum`
    For more advanced cases, use Fitroom : https://github.com/radis/fitroom
    """

    # Get initial values of fitted parameters
    model_input = fit_parameters

    compute_los_model = lambda model_input: model(factory, model_input)

    # Calculate initial Spectrum, by showing all steps.
    # old_verbose = factory.verbose
    factory.verbose = 0  # reduce verbose during calculation.
    s0 = compute_los_model(model_input)  # Load enrgies; initialize all caches, etc.
    factory.verbose = 3  # increase verbose level for more details.
    s0 = compute_los_model(
        model_input
    )  # .plot('transmittance_noslit', nfig='Initial spectrum')
    factory.verbose = 0  # reduce verbose during calculation.

    # %% Leastsq version

    # %%
    # User Params
    # -----------

    fit_params = list(fit_parameters.keys())
    bounds_arr = np.array([bounds[k] for k in fit_params])
    # fit_units = ['K', 'K', 'K']
    import astropy.units as u

    fit_units = []
    for k in fit_params:
        if isinstance(fit_parameters[k], u.Quantity) or isinstance(
            fit_parameters[k], u.Unit
        ):
            fit_units.append(u.Unit(fit_parameters[k]))
        else:
            fit_units.append("")

    if len(s_exp.get_vars()) != 1:
        raise ValueError(
            "More than one spectral array in experimental spectrum"
            + f"({len(s_exp.get_vars())}) : {s_exp.get_vars()}. "
            + "Choose one only with `s_exp.take()`"
        )
    fit_variable = s_exp.get_vars()[0]

    # -----------------------------------------------------------------------------
    # -----------------------------------------------------------------------------
    #                 FITTING MACHINERY    (you shouldnt need to edit this)
    #                  ... just a few functions to make nice plots along
    #                  ... the fitting procedure
    # -----------------------------------------------------------------------------
    # -----------------------------------------------------------------------------

    # %%
    # Algorithm Params
    # ----------------

    history_x = []
    history_res = []
    maxiter = maxiter

    def print_fit_values(fit_values):
        return ",".join(
            [
                "{0}={1}{2}".format(
                    fit_params[i], np.round(fit_values[i], 0), fit_units[i]
                )
                for i in range(len(fit_params))
            ]
        )

    def generate_spectrum(fit_values):

        # Generate dictionary
        inputs = model_input.copy()
        for k, v in zip(fit_params, fit_values):
            inputs[k] = v

        # Calculate the theoretical model
        s = compute_los_model(inputs)

        return s

    if plot:
        # Graph with plot diff
        # figSpec, axSpec = plt.subplots(num='diffspectra')
        figSpec, axSpec = plot_diff(s_exp, s0, fit_variable, nfig="diffspectra")
        lineSpec = axSpec[0].get_lines()[1]
        lineDiff = axSpec[1].get_lines()[0]

    def cost_function(fit_values, plot=None):
        """ Return error on Spectrum s vs experimental spectrum"""

        s = generate_spectrum(fit_values)

        # Delete unecessary variables (for a faster resampling)
        for var in [k for k in s._q.keys() if k not in [fit_variable, "wavespace"]]:
            del s._q[var]

        if plot:
            s_diff = get_diff(s_exp, s, var=fit_variable)
            lineSpec.set_data(s.get(fit_variable))
            lineDiff.set_data(s_diff)
            figSpec.canvas.draw()
            # plt.pause(0.001)
            plt.pause(0.05)

        s.resample(s_exp, energy_threshold=2e-2)

        return get_residual(s, s_exp, fit_variable, ignore_nan=True, norm="L2")

    def log_cost_function(fit_values, plot=None):
        """ Calls the cost_function, and write the values to the Log history  """

        res = cost_function(fit_values, plot=plot)

        history_x.append(fit_values)
        history_res.append(res)

        return res

    # Graph with residual
    # ... unlike 1D we cant plot the temperature here. Just plot the iteration

    if plot:
        plt.close("residual")
        figRes, axRes = plt.subplots(num="residual", figsize=(13.25, 6))
        axValues = axRes.twinx()
        # ax = fig.gca()
        # ax.set_ylim((bounds[0]))

    fit_values_min, fit_values_max = bounds_arr.T
    res0 = log_cost_function(fit_values_min, plot=plot)
    res1 = log_cost_function(fit_values_max, plot=plot)

    if plot:
        (lineRes,) = axRes.plot((1, 2), (res0, res1), "-ko")
        (lineLast,) = axRes.plot(2, res0, "or")  # last iteration in red
        lineValues = {}
        for i, k in enumerate(fit_params):
            lineValues[k] = axValues.plot(
                (1, 2), (fit_values_min[i], fit_values_max[i]), "-", label=k
            )[0]
        axRes.set_xlim((0, maxiter))
        axRes.set_ylim(ymin=0)
        axRes.set_xlabel("Iteration")
        axRes.set_ylabel("Residual")
        figRes.legend()

    factory.verbose = False
    factory.warnings["NegativeEnergiesWarning"] = "ignore"

    global ite
    ite = 2

    if plot:
        plot_every = max(
            1, int(0.2 / s0.conditions["calculation_time"])
        )  # refresh plot every X calculations

    def cost_and_plot_function(fit_values):
        """Return error on Spectrum s vs experimental spectrum

        This is the function that is called by minimize()"""
        global ite
        ite += 1
        # Plot one spectrum every X ites
        plot_ite = plot and not ite % plot_every

        res = log_cost_function(fit_values, plot=plot_ite)

        if plot:
            # Add to plot history
            x, y = lineRes.get_data()
            lineRes.set_data((np.hstack((x, ite)), np.hstack((y, res))))
            # Add values to history
            for k, v in zip(fit_params, fit_values):
                x, y = lineValues[k].get_data()
                lineValues[k].set_data((np.hstack((x, ite)), np.hstack((y, v))))
            # Plot last
            lineLast.set_data((ite, res))

            figRes.canvas.draw()

        print(
            "{0}, Residual: {1:.4f} {2}".format(
                print_fit_values(fit_values),
                res,
                " 🏆" if res == min(history_res) else "",
            ),
            flush=True,
        )

        return res

    # >>> This is where the fitting loop happens
    print("\nNow starting the fitting process:")
    print("---------------------------------\n")
    best = minimize(
        cost_and_plot_function,
        (fit_values_max + fit_values_min) / 2,
        # method='L-BFGS-B',
        method="TNC",
        jac=None,
        bounds=bounds_arr,
        options={
            "maxiter": int(max(1, maxiter / 4)),  # somehow.
            "eps": 20,
            #                         'ftol':1e-10,
            # 'gtol':1e-10,
            "disp": True,
        },
    )
    # <<<

    s_best = generate_spectrum(best.x)

    if best.success:
        print("Final {0}: {1}{2}".format(fit_params, np.round(best.x), fit_units))

    # Res history

    # ... what does history say:
    print(
        "Best: {0}: {1}{2} reached at iteration {3}/{4}".format(
            fit_params,
            history_x[np.argmin(history_res)],
            fit_units,
            np.argmin(history_res),
            best.nfev,
        )
    )

    # ... note that there are more function evaluations (best.nfev) that actual solver
    # ... iterations (best.nit) because the Jacobian is calculated numerically with
    # ... internal function calls

    return s_best


if __name__ == "__main__":

    # %% Get Fitted Data

    # Data from Dang, adapted by Klarenaar, digitized by us
    s_exp = Spectrum.from_txt(
        getValidationCase(
            join(
                "test_CO2_3Tvib_vs_klarenaar_data", "klarenaar_2017_digitized_data.csv"
            )
        ),
        "transmittance_noslit",
        waveunit="cm-1",
        unit="",
        delimiter=",",
        name="Klarenaar 2017",
    )

    # %% Calculate

    sf = SpectrumFactory(
        2284.2,
        2284.6,
        wstep=0.001,  # cm-1
        pressure=20 * 1e-3,  # bar
        db_use_cached=True,
        lvl_use_cached=True,
        cutoff=1e-25,
        isotope="1,2",
        path_length=10,  # cm-1
        mole_fraction=0.1 * 28.97 / 44.07,
        broadening_max_width=1,  # cm-1
        medium="vacuum",
        export_populations=None,  # 'vib',
    )
    sf.warnings["MissingSelfBroadeningWarning"] = "ignore"
    setup_test_line_databases()
    sf.load_databank("HITEMP-CO2-TEST")

    s_best = sf.fit_spectrum(
        s_exp.take("transmittance_noslit"),
        model=Tvib12TrotModel,
        fit_parameters={
            "T12": 517,
            "T3": 2641,
            "Trot": 491,
        },
        bounds={"T12": [300, 2000], "T3": [300, 5000], "Trot": [300, 2000]},
        plot=True,
        maxiter=120,
    )
    plot_diff(s_exp, s_best)
