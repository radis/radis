# -*- coding: utf-8 -*-
"""
Functions and constants used in :class:`~radis.spectrum.spectrum.Spectrum`
object

-------------------------------------------------------------------------------


"""

import subprocess

import numpy as np

from radis.misc.basics import partition

# %% Definitions

# Waverange, wavespaces
# ... aliases:
WAVENUM_UNITS = ["cm", "cm-1", "wavenumber"]
WAVELEN_UNITS = ["nm", "wavelength"]
WAVELENVAC_UNITS = ["nm_vac", "nm_vacuum"]
# ... internal names
WAVESPACE = ["nm", "nm_vac", "cm-1"]
""" list: wavespace:
- ``'nm'``: wavelength (in air)
- ``'nm_vac'``: wavelength (in vacuum)
- ``'cm-1'``: wavenumber"""

# Spectral quantities
CONVOLUTED_QUANTITIES = ["radiance", "transmittance", "emissivity"]
"""list: name of spectral quantities after convolution with slit function"""
NON_CONVOLUTED_QUANTITIES = [
    "radiance_noslit",
    "transmittance_noslit",
    "emisscoeff",
    "emisscoeff_continuum",
    "absorbance",
    "abscoeff",  # opacity
    "abscoeff_continuum",
    "emissivity_noslit",
]
"""list: name of spectral quantities not convolved with slit function

See :ref:`the description of spectral quantities <label_spectral_quantities>`"""
SPECTRAL_QUANTITIES = CONVOLUTED_QUANTITIES + NON_CONVOLUTED_QUANTITIES
"""list: all spectral quantities defined in a :class:`~radis.spectrum.spectrum.Spectrum`
object.

See :ref:`the description of spectral quantities <label_spectral_quantities>`"""

# note: it is hardcoded (and needed) that quantities that are convoluted are
# generated from a non convoluted quantity with the same name + _noslit
for _q in CONVOLUTED_QUANTITIES:
    assert _q + "_noslit" in NON_CONVOLUTED_QUANTITIES

# %% Special names when printing conditions

# These parameters are shown below "Physical Conditions" rather than "Computation
# Parameters
PHYSICAL_PARAMS = [
    "molecule",
    "wavenum_max",
    "wavenum_min",
    "mole_fraction",
    "isotope",
    "state",
    "path_length",
    "medium",
    "self_absorption",
    "slit_function_base",
    "pressure_mbar",
    "wavelength_min",
    "wavelength_max",
    "Telec",
    "Tvib",
    "Trot",
    "Tgas",
    "vib_distribution",
    "rot_distribution",
    "overpopulation",
    "thermal_equilibrium",
]
"""list: physical conditions under which the Spectrum was calculated/measured.
When printing an object, these parameters are shown below "Physical Conditions"
If a parameter is not in this list, it is either in "Computation Parameters"
(non-physical parameters that can have an influence on the Spectrum, e.g, cutoffs
and thresholds) or in "Informative Params" (descriptive parameters that have absolutely no
impact on the spectrum, e.g, number of lines calculated or calculation time)
"""

# Informative parameters
# ... Parameters that should be saved in the Spectrum objects, but ignored when
# ... comparing two spectra.
# ... should be written here only these parameters that cannot affect the
# ... physical result. In particular, all parameters relative to performance
# ... should be added here.
INFORMATIVE_PARAMS = [
    "db_use_cached",
    "lvl_use_cached",
    "chunksize",
    "calculation_time",
    "total_lines",
    "lines_calculated",
    "lines_cutoff",
    "lines_in_continuum",
    "Nprocs",
    "warning_broadening_threshold",
    "warning_linestrength_cutoff",
    "load_energies",
    "export_lines",
    "export_populations",
]
""" list: Informative parameters. Parameters that should be saved in the Spectrum
objects, but ignored when comparing two spectra. Should be written here only
these parameters that cannot affect the physical result. In particular, all
parameters relative to performance should be added here.

Notes
-----

units for these parameters are stored in Spectrum.cond_units and are defined
by the generating class (ex: SpectrumFactory)"""

# %% Util functions


def cast_waveunit(unit, force_match=True):
    """Standardize unit formats."""
    if unit in WAVELEN_UNITS:
        return "nm"
    if unit in WAVELENVAC_UNITS:
        return "nm_vac"
    elif unit in WAVENUM_UNITS:
        return "cm-1"
    elif force_match:
        raise ValueError(
            "Unknown wavespace unit: {0}. Should be one of {1}".format(
                unit, WAVELEN_UNITS + WAVELENVAC_UNITS + WAVENUM_UNITS
            )
        )
    else:
        return unit  # dont convert


def make_up(label):
    """Cosmetic changes on label, before plot.

    Parameters
    ----------

    label: str
    """

    # Improve units
    label = label.replace(r"cm-1", r"cm⁻¹")
    label = label.replace(r"m^-1", r"m⁻¹")
    label = label.replace(r"m2", r"m²")
    label = label.replace(r"m3", r"m³")
    label = label.replace(r"I/I0", r"I/I₀")  # transmittance unit

    # Improve text
    if not "_noslit" in label:
        # small hack: capitalize, to make a difference with non slit value
        label = label.replace(
            "radiance", "Radiance"
        )  # make a small difference between no_slit and slit while plotting
        label = label.replace(
            "transmittance", "Transmittance"
        )  # make a small difference between no_slit and slit while plotting
    # ... Remove _noslit
    label = label.replace("_noslit", "")

    return label


def make_up_unit(Iunit, var):
    """Additional cosmetic changes for units on label, before plot.

    Parameters
    ----------

    Iunit: str
        input unit

    var: str
        spectral variable. Ex: ``transmittance``
    """
    Iunit = Iunit.replace(r"um", r"µm")
    Iunit = Iunit.replace(r"cm-1", r"cm⁻¹")  # welcome to unicode ! ;)
    Iunit = Iunit.replace(r"m2", r"m²")

    if Iunit == "":
        # give more explicit unit for the user:
        if var in ["transmittance", "transmittance_noslit"]:
            Iunit = r"I/I0"
        elif var == "absorbance":
            Iunit = r"-ln(I/I0)"
        elif var in ["emissivity_no_slit", "emissivity"]:
            Iunit = r"ε"
        elif var in ["radiance", "radiance_noslit"]:
            Iunit = r"norm"

    return Iunit


def format_xlabel(wunit, plot_medium):
    """Used by :py:meth:`radis.spectrum.spectrum.Spectrum.plot` and
    :py:func:`radis.spectrum.compare.plot_diff`

    Parameters
    ----------

    wunit: ``'default'``, ``'nm'``, ``'cm-1'``, ``'nm_vac'``,
        wavelength air, wavenumber, or wavelength vacuum. If ``'default'``,
        Spectrum :py:meth:`~radis.spectrum.spectrum.Spectrum.get_waveunit` is used.

    plot_medium: bool, ``'vacuum_only'``
        if ``True`` and ``wunit`` are wavelengths, plot the propagation medium
        in the xaxis label (``[air]`` or ``[vacuum]``). If ``'vacuum_only'``,
        plot only if ``wunit=='nm_vac'``. Default ``'vacuum_only'``
        (prevents from inadvertently plotting spectra with different propagation
        medium on the same graph).

    """
    if wunit == "cm-1":
        xlabel = "Wavenumber (cm-1)"
    elif wunit == "nm":
        if plot_medium and plot_medium != "vacuum_only":
            xlabel = "Wavelength [air] (nm)"
        else:
            xlabel = "Wavelength (nm)"
    elif wunit == "nm_vac":
        if plot_medium and plot_medium != "vacuum_only":
            xlabel = "Wavelength [vacuum] (nm)"
        else:
            xlabel = "Wavelength (nm)"
    else:
        raise ValueError(wunit)

    return make_up(xlabel)


def print_conditions(
    conditions,
    units,
    phys_param_list=PHYSICAL_PARAMS,
    info_param_list=INFORMATIVE_PARAMS,
):
    """Print all Spectrum calculation parameters.

    Parameters
    ----------

    phys_param_list: list
        These parameters are shown below "Physical Conditions" rather than "Computation
        Parameters. See :data:`~radis.spectrum.utils.PHYSICAL_PARAMS` for more
        information.

    info_param_list: list
        These parameters are shown below "Information" rather than "Computation
        Parameters. See :data:`~radis.spectrum.utils.INFORMATIVE_PARAMS` for more
        information.

    See Also
    --------

    :data:`~radis.spectrum.utils.PHYSICAL_PARAMS`, :data:`~radis.spectrum.utils.INFORMATIVE_PARAMS`
    """

    def align(a, space=20):
        """fix alignement."""
        return a + " " * max(1, (space - len(str(a))))

    def print_param(k):
        """Special formatting for nicely printing conditions."""
        v_k = conditions[k]
        # Add extra arguments based on arbitrary conditions
        args = []
        if k in units:
            args.append(units[k])
        # ... fill here for other args

        # Special formatting
        try:
            if k in [
                "wavenum_max_calc",
                "wavenum_min_calc",
                "wavelength_max",
                "wavelength_min",
                "wavenum_max",
                "wavenum_min",
            ]:
                v_k_str = "{0:.4f}".format(v_k)
            elif k in ["lines_calculated", "lines_in_continuum", "lines_cutoff"]:
                # Add comma separator for thousands
                v_k_str = "{0:,d}".format(v_k)
            else:
                # Default to printing str
                v_k_str = "{0}".format(v_k)
        except ValueError:
            # Default to printing str
            v_k_str = "{0}".format(v_k)

        # Crop
        if len(v_k_str) > 102:  # cut if too long
            v_k_str = v_k_str[:100] + "..."
        print(" " * 2, align(k), v_k_str, *args)

    phys_param, non_phys_param = partition(lambda x: x in phys_param_list, conditions)

    info_param, non_phys_param = partition(
        lambda x: x in info_param_list, non_phys_param
    )

    print("Physical Conditions")
    print("-" * 40)
    for k in sorted(phys_param):
        print_param(k)

    print("Computation Parameters")
    print("-" * 40)
    for k in sorted(non_phys_param):
        print_param(k)

    if len(info_param) > 0:
        print("Information")
        print("-" * 40)
        for k in sorted(info_param):
            print_param(k)

    print("-" * 40)

    # print gas_inp (information on each gas slab) if exists (specifically
    # for Specair output)
    if "gas_inp" in conditions:
        try:
            for slab in conditions["gas_inp"]:
                print("Slab", slab)
                slab.print_conditions()
        except:
            pass

    return None


def has_nan(s):
    """

    Parameters
    ----------
    s : Spectrum
        radis Spectrum.

    Returns
    -------
    b : bool
        returns whether Spectrum has ``nan``

    Note
    ----

    ``print(s)`` will also show which spectral quantities have ````nan.

    """

    for k, v in s._get_items().items():
        if np.isnan(v[1]).any():
            return True
    return False


# %% Plot helper


def split_and_plot_by_parts(w, I, *args, **kwargs):
    from warnings import warn

    from radis.misc.plot import split_and_plot_by_parts

    warn(
        DeprecationWarning(
            "split_and_plot_by_parts() moved in radis.misc.plot in 0.9.30"
        )
    )
    return split_and_plot_by_parts(w, I, *args, **kwargs)


# %% Profiler

#%%
def dict_to_tree(pro, name):
    """Convert to ::
    Tree =
    {"name": str,
     "value":float,
     "children":[list of Tree]}

    """
    if isinstance(pro, dict):
        new_dict = {"value": pro["value"], "name": name}
        children = {k: v for k, v in pro.items() if k != "value"}
        if len(children) > 0:
            new_dict["children"] = [
                dict_to_tree(v, name=k) for k, v in children.items()
            ]
        return new_dict
    else:
        return {"value": pro, "name": name}

    return new_dict


#%%


def generate_perf_profile(s):
    """Visual/interactive performance profile

    See typical output in https://github.com/radis/radis/pull/324

    Parameters
    ----------
    s: Spectrum

    See Also
    --------
    :py:meth:`~radis.spectrum.generate_perf_profile`
    """
    from pstats import Stats

    st = Stats()

    """
    From Lib/profile.py:Profile

        Timing data for each function is stored as a 5-tuple in the dictionary
        self.timings[].  The index is always the name stored in self.cur[-3].
        The following are the definitions of the members:

        [0] = The number of times this function was called, not counting direct
              or indirect recursion,
        [1] = Number of times this function appears on the stack, minus one
        [2] = Total time spent internal to this function
        [3] = Cumulative time that this function was present on the stack.  In
              non-recursive functions, this is the total execution time from start
              to finish of each invocation of a function, including time spent in
              all subfunctions.
        [4] = A dictionary indicating for each function name, the number of times
              it was called by us.
    """

    def parse_profiler(tree, parent=None, parent_time=None):
        """

        Parameters
        ----------
        tree : dict
            ::
            {"name": str,
             "value":float,
             "children":[list of Tree]}
        parent : tuple (str, int, str)
            filename(str), line_number(int), parent function name(str)
        parent_time : int, int,
            ``cc, nc, tt, ct`` of parent function.
            - The number of times this function was called, not counting direct
              or indirect recursion,
            - Number of times this function appears on the stack, minus one
            - Total time spent internal to this function
            - Cumulative time that this function was present on the stack.

        """
        func_name = (
            "",
            1,
            tree["name"],
        )  # normally : filename(str), line_number(int), function name (str)

        # Parse children
        if "children" in tree:
            children_cumtime = sum([child["value"] for child in tree["children"]])
            for child in tree["children"]:
                # child_name = "", 1, child["name"]
                child_time = child["value"]
                parse_profiler(
                    child, parent=func_name, parent_time=(1, 0, child_time, child_time)
                )

            if parent:
                parent_time = (
                    parent_time[0],
                    parent_time[1],
                    parent_time[2] - children_cumtime,
                    parent_time[3],
                )
            tt = tree["value"] - children_cumtime
            ct = tree["value"]

        else:
            tt = tree["value"]
            ct = tree["value"]

        # Store value
        st.stats[func_name] = (1, 0, tt, ct, {parent: parent_time} if parent else {})

    profiler = s.conditions["profiler"].copy()
    # Add total calculation time:
    profiler.update({"value": s.conditions["calculation_time"]})

    # Fix
    if "spectrum_calc_before_obj" in profiler:
        from warnings import warn

        warn("`spectrum_calc_before_obj` still in profiler keys")
        del profiler["spectrum_calc_before_obj"]

    perf_tree = dict_to_tree(profiler, name="calculation_time")

    parse_profiler(perf_tree)
    st.dump_stats("spectrum.prof")

    subprocess.Popen(["tuna", "spectrum.prof"])
