# -*- coding: utf-8 -*-
"""
Functions and constants used in :class:`~radis.spectrum.spectrum.Spectrum`
object

-------------------------------------------------------------------------------


"""

import matplotlib.pyplot as plt
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
    """Plot two discontinued arrays (typically a spectrum) without showing
    junctions: first identify junctions then split and plot separately.

    Useful for plotting an experimental spectrum defined on different, non overlapping
    ranges without showing connecting lines between the ranges, or to plot an
    experimental spectrum defined on overlapping ranges, without showing connecting
    lines neither.

    Parameters
    ----------

    w, I: arrays
        typically output of :py:func:`~numpy.hstack`.

    Other Parameters
    ----------------

    split_threshold: int
        number of standard deviation for threshold. Default 10

    ax: matplotlib axe
        plot on a particular axe

    kwargs: dict
        forwarded to :func:`~matplotlib.pyplot.plot`

    cutwings: int
        discard elements on the side. Default 0
    """

    from publib.tools import keep_color

    # Get defaults
    ax = kwargs.pop("ax", None)
    if "ax" == None:
        ax = plt.gca()
    split_threshold = kwargs.pop("split_threshold", 10)  # type: int
    cutwings = kwargs.pop("cutwings", 0)  # type: int
    label = kwargs.pop("label", None)  # type: str

    # identify joints
    dw = np.diff(w)
    dwmean = dw.mean()
    joints = np.argwhere((abs(dw - dwmean) > split_threshold * dw.std())) + 1

    # Split
    if len(joints) > 0:
        ws = np.split(w, joints.flatten())
        Is = np.split(I, joints.flatten())

        # Plot separately
        out = []
        for i, (wi, Ii) in enumerate(zip(ws, Is)):
            if cutwings:
                wi = wi[cutwings:-cutwings]
                Ii = Ii[cutwings:-cutwings]
            if i == 0:  # label once only
                out.append(
                    ax.plot(
                        wi, Ii, *args, **dict(list(kwargs.items()) + [("label", label)])
                    )
                )
            else:
                keep_color()
                out.append(ax.plot(wi, Ii, *args, **kwargs))

        return list(zip(*out))

    else:
        if cutwings:
            w = w[cutwings:-cutwings]
            I = I[cutwings:-cutwings]
        return ax.plot(w, I, *args, **dict(list(kwargs.items()) + [("label", label)]))
