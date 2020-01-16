# -*- coding: utf-8 -*-
"""
Functions and constants used in :class:`~radis.spectrum.spectrum.Spectrum` 
object

-------------------------------------------------------------------------------


"""

from __future__ import print_function, absolute_import, division, unicode_literals
from radis.misc.basics import partition

# %% Definitions

# Waverange, wavespaces
# ... aliases:
WAVENUM_UNITS = ["cm", "cm-1", "cm_1", "wavenumber"]
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
    "abscoeff",
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
    "parallel",
    "Nprocs",
    "Ngroups",
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
    """ Standardize unit formats """
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
    """ Cosmetic changes on label, before plot 


    Parameters    
    ----------

    label: str
    """

    # Improve units
    label = label.replace(r"cm_1", r"cm-1")
    label = label.replace(r"cm-1", r"cm$^\mathregular{-1}$")
    label = label.replace(r"m^-1", r"m$^\mathregular{-1}$")
    label = label.replace(r"m2", r"m$^\mathregular{2}$")
    label = label.replace(r"m3", r"m$^\mathregular{3}$")
    label = label.replace(r"I/I0", r"I/I$_\mathregular{0}$")  # transmittance unit

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


def format_xlabel(wunit, plot_medium):
    """ Used by :py:meth:`radis.spectrum.spectrum.Spectrum.plot` and 
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
    """ Print all Spectrum calculation parameters 

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
        """ fix alignement """
        return a + " " * max(1, (space - len(str(a))))

    def print_param(k):
        """ Special formatting for nicely printing conditions """
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
