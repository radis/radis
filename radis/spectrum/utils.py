# -*- coding: utf-8 -*-
"""
Functions and constants used in :class:`~radis.spectrum.spectrum.Spectrum`
object

-------------------------------------------------------------------------------


"""

import subprocess

from radis.misc.basics import partition

# %% Definitions

# Waverange, wavespaces
# ... aliases:
WAVENUM_UNITS = ["cm", "cm-1"]
WAVELEN_UNITS = ["nm", "nm_air"]
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
    "xsection",  # cross-sections
]
"""list: name of spectral quantities not convolved with slit function

See :ref:`the description of spectral arrays <label_spectral_arrays>`"""
SPECTRAL_QUANTITIES = CONVOLUTED_QUANTITIES + NON_CONVOLUTED_QUANTITIES
"""list: all spectral quantities defined in a :class:`~radis.spectrum.spectrum.Spectrum`
object.

See :ref:`the description of spectral arrays <label_spectral_arrays>`"""

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
    # parameters added to reduce default verbose :
    "levelsfmt",
    "wavenum_max_calc",
    "wavenum_min_calc",
    "dxG",  # TODO : keep it a Computational parameter, but of verbose-level 2 so it's not shown by default
    "dxL",
    "export_rovib_fraction",
    "parfuncpath",
]
""" list: Informative parameters. Parameters that should be saved in the Spectrum
objects, but ignored when comparing two spectra. Should be written here only
these parameters that cannot affect the physical result. In particular, all
parameters relative to performance should be added here.

Notes
-----

units for these parameters are stored in Spectrum.cond_units and are defined
by the generating class (ex: SpectrumFactory)"""

CONFIG_PARAMS = [
    "GRIDPOINTS_PER_LINEWIDTH_WARN_THRESHOLD",
    "GRIDPOINTS_PER_LINEWIDTH_ERROR_THRESHOLD",
    "SPARSE_WAVERANGE",
    "DEFAULT_DOWNLOAD_PATH",
]
""" list: these parameters are read from radis.config and stored in the Spectrum
objects. Should be added here only the parameters that may have an impact on the computation,
for instance the one defining the 'auto' thresholds
"""

# %% Util functions


def cast_waveunit(unit, force_match=True):
    """Standardize unit formats, return either "nm", "nm_vac" or "cm-1"."""
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

    # Hardfix for missing characters
    SUPERSCRIPT_MINUS = "⁻"
    from matplotlib.font_manager import FontProperties, findfont

    if "arial" in findfont(FontProperties(family=["sans-serif"])):
        SUPERSCRIPT_MINUS = "$^{-}$"

    # Improve units
    label = label.replace(r"cm-1", r"cm" + f"{SUPERSCRIPT_MINUS}" + r"¹")
    label = label.replace(r"m^-1", r"m" + f"{SUPERSCRIPT_MINUS}" + r"¹")
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

    label = label.replace("xsection", "cross-section")

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
    config_param_list=CONFIG_PARAMS,
    verbose=2,
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
    config_param_list: list
        These parameters are read from radis.config file. See
        :data:`~radis.spectrum.utils.CONFIG_PARAMS` for more
        information.
    verbose: int
        if ``1`` or ``True``, only physical and computational parameters are shown. If ``2``,
        all parameters (including config & informative) are shown. Default ``2``.

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

    config_param, non_phys_param = partition(
        lambda x: x in config_param_list, non_phys_param
    )

    # TODO : add verbose level so that Trot, Tvib are not shown if equilibrium
    # and verbose<2 . Same for rot_distriubtion; vib_distriubtion

    print("Physical Conditions")
    print("-" * 40)
    for k in sorted(phys_param):
        print_param(k)

    print("Computation Parameters")
    print("-" * 40)
    for k in sorted(non_phys_param):
        print_param(k)

    if verbose >= 2:
        print("Config parameters")
        print("-" * 40)
        for k in sorted(config_param):
            print_param(k)

        if len(info_param) > 0:
            print("Information")
            print("-" * 40)
            for k in sorted(info_param):
                print_param(k)

    print("-" * 40)

    if verbose >= 2:

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
    """
    Parameters
    ----------
    pro: dict
        of the form::
            {"value":float,    # keyword 'value' is expected
             "some_key":float
             "some_key2":float,
             "some_key3":dict  # nested dict of the same form
             "some_key4":dict}
    name: str

    Returns
    -------
    dict: of the form ::
        Tree =
        {"name": str,
         "value":float,
         "children":[list of Tree]}

    See Also
    --------
    Used in :py:func:`~radis.spectrum.utils.print_perf_profile` and
    :py:func:`~radis.spectrum.utils.generate_perf_profile`

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


#%%


def print_perf_profile(
    profiler,
    total_time=None,
    number_format="{:.3f}",
    precision=16,
    first_line="profiler :",
):
    r"""Prints Profiler output dictionary in a structured manner.

    Parameters
    ----------
    profiler: dict
        of the form::
            {"value":float,    # keyword 'value' is expected
             "some_key":float
             "some_key2":float,
             "some_key3":dict  # nested dict of the same form
             "some_key4":dict}

    Other Parameters
    ----------------
    total_time: float
        total calculation time. If ``None``, take the key ``"value"`` of ``profiler``
    precision: int, optional
        total number of blocks. Default 16.

    Example
    -------
    ::

        Spectrum.print_perf_profile()

        # output >>
            spectrum_calculation      0.189s ████████████████
                check_line_databank              0.000s
                check_non_eq_param               0.042s ███
                fetch_energy_5                   0.015s █
                calc_weight_trans                0.008s
                reinitialize                     0.002s
                    copy_database                    0.000s
                    memory_usage_warning             0.002s
                    reset_population                 0.000s
                calc_noneq_population            0.041s ███
                    part_function                    0.035s ██
                    population                       0.006s
                scaled_non_eq_linestrength       0.005s
                    map_part_func                    0.001s
                    corrected_population_se          0.003s
                calc_emission_integral           0.006s
                applied_linestrength_cutoff      0.002s
                calc_lineshift                   0.001s
                calc_hwhm                        0.007s
                generate_wavenumber_arrays       0.001s
                calc_line_broadening             0.074s ██████
                    precompute_LDM_lineshapes        0.012s
                    LDM_Initialized_vectors          0.000s
                    LDM_closest_matching_line        0.001s
                    LDM_Distribute_lines             0.001s
                    LDM_convolve                     0.060s █████
                    others                           0.001s
                calc_other_spectral_quan         0.003s
                generate_spectrum_obj            0.000s
                others                           -0.016s

    See Also
    --------
    :py:meth:`~radis.spectrum.spectrum.Spectrum.print_perf_profile`,
    :py:meth:`~radis.lbl.factory.SpectrumFactory.print_perf_profile`
    """
    PRECISION = precision  #  number of blocks █
    TAB = 4  # number of spaces for indentation

    if total_time is None:
        total_time = profiler["value"]

    def scale(time):
        # return '|'*int(time/total_time*10)
        return "\u2588" * int(time / total_time * PRECISION)
        # '\u2588'  is  █

    def walk_print_tree(prof, name, level, write_number_column):
        if "value" in prof:
            text = " " * TAB * level + name
            fill_spaces = " " * (write_number_column - len(text))
            print(
                text,
                fill_spaces,
                "" + number_format.format(prof["value"]) + "s",
                scale(prof["value"]),
            )
        max_name_length = max(len(k) for k in prof.keys() if k != "value")
        write_number_column = max(
            write_number_column, len(" " * TAB * (level + 1)) + max_name_length
        )
        write_number_column += TAB

        total_time = 0
        for k, v in prof.items():
            if k == "value":
                pass
            elif isinstance(v, dict):
                total_time += walk_print_tree(
                    v,
                    name=k,
                    level=level + 1,
                    write_number_column=write_number_column,
                )
            elif isinstance(v, float):
                text = " " * TAB * (level + 1) + k
                fill_spaces = " " * (write_number_column - len(text))
                print(text, fill_spaces, "" + number_format.format(v) + "s", scale(v))
                total_time += v
            else:
                raise ValueError(type(v))

        # print missing time / self-time
        if "value" in prof:
            missing_time = prof["value"] - total_time
            if float(number_format.format(missing_time)) != 0:
                # we dont add 0 numbers
                text = " " * TAB * (level + 1) + "others"
                fill_spaces = " " * (write_number_column - len(text))
                print(
                    text,
                    fill_spaces,
                    "" + number_format.format(missing_time) + "s",
                    scale(missing_time),
                )

        return total_time

    print(first_line)
    walk_print_tree(profiler, name="", level=0, write_number_column=0)


def generate_perf_profile(profiler):
    """Visual/interactive performance profile

    Requires ``tuna`` to be installed.

    See typical output in https://github.com/radis/radis/pull/325

    .. image:: https://user-images.githubusercontent.com/16088743/128018032-6049be72-1881-46ac-9d7c-1ed89f9c4f42.png
        :alt: https://user-images.githubusercontent.com/16088743/128018032-6049be72-1881-46ac-9d7c-1ed89f9c4f42.png
        :target: https://user-images.githubusercontent.com/16088743/128018032-6049be72-1881-46ac-9d7c-1ed89f9c4f42.png


    .. note::
        You can also profile with `tuna` directly::

            python -m cProfile -o program.prof your_radis_script.py
            tuna your_radis_script.py

    Parameters
    ----------
    profiler: dict
        of the form::
            {"value":float,    # keyword 'value' is expected
             "some_key":float
             "some_key2":float,
             "some_key3":dict  # nested dict of the same form
             "some_key4":dict}

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

    perf_tree = dict_to_tree(profiler, name="calculation_time")

    parse_profiler(perf_tree)
    st.dump_stats("spectrum.prof")

    return subprocess.Popen(["tuna", "spectrum.prof"])


if __name__ == "__main__":
    from radis.test.spectrum.test_utils import test_perf_profile

    test_perf_profile()
