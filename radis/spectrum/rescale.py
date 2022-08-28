# -*- coding: utf-8 -*-
"""Functions to update :class:`~radis.spectrum.spectrum.Spectrum` objects with
new spectral quantities that can be derived from existing ones, or rescale
path_length or mole_fraction, or add overpopulations.

Most of these are binded as methods to the Spectrum class, but stored here to
unload the spectrum.py file


Equations derived from the code using `pytexit <https://pytexit.readthedocs.io/en/latest/>`__

-------------------------------------------------------------------------------
"""

from warnings import warn

import astropy.units as u
import numpy as np
from numpy import exp
from numpy import log as ln

from radis.misc.basics import all_in, any_in, compare_lists
from radis.misc.debug import printdbg
from radis.phys.units_astropy import convert_and_strip_units
from radis.spectrum.equations import abscoeff2xsection, calc_radiance
from radis.spectrum.utils import CONVOLUTED_QUANTITIES, NON_CONVOLUTED_QUANTITIES

# List of all spectral variables sorted by priority during recomputation
# (ex: first get abscoeff, then try to calculate emisscoeff, etc.)
ordered_keys = [
    "abscoeff",
    "emisscoeff",
    "absorbance",
    "radiance_noslit",
    "transmittance_noslit",
    "emissivity",
    "emissivity_noslit",
    "radiance",
    "transmittance",
    "xsection",
]
"""list: List of all spectral variables sorted by priority during recomputation
(ex: first get abscoeff, then try to calculate emisscoeff, etc.)
"""

# ... variables that cannot be rescaled (or not implemented):
non_rescalable_keys = ["abscoeff_continuum", "emisscoeff_continuum"]
"""str:  variables that cannot be rescaled (or not implemented): """
# ... Check we have everyone (safety check!):
# ... if it fails here, then we may have added a new key without adding a scaling
# ... method. Explicitely add it in non_rescalableçkeys so an error is raised
# ... if trying to rescale a Spectrum that has such a quantity
assert (
    compare_lists(
        ordered_keys + non_rescalable_keys,
        CONVOLUTED_QUANTITIES + NON_CONVOLUTED_QUANTITIES,
        verbose=False,
    )
    == 1
)


def _build_update_graph(
    spec,
    optically_thin=None,
    equilibrium=None,
    path_length=None,
    p_T_x=None,
    no_change=False,
):
    """Find inheritances properties (dependencies and equivalences) between all
    spectral variables based on the spectrum conditions (equilibrium, optically
    thin, known path length?)

    Parameters
    ----------
    spec: Spectrum
        a :class:`~radis.spectrum.spectrum.Spectrum` object

    Other Parameters
    ----------------
    optically_thin: boolean
        know whether the Spectrum should be considered optically thin to build
        the equivalence graph tree. If None, the value stored in the Spectrum is used.
        Default ``None``
    equilibrium: boolean
        know whether the Spectrum should be considered at equilibrium to build
        the equivalence graph tree. If None, the value stored in the Spectrum is used.
        Default ``None``
    path_length: boolean
        know whether the path length is given to build the equivalence graph tree.
        If None, ``path_length`` is looked up in the Spectrum condition. Default ``None``
    p_T_x: boolean
        whether all of the temperature, pressure and mole fractions are known. If None,
        they are looked up in the Spectrum conditions. Default ``None``
    no_change: boolean
        if True, signals that we are somehow rescaling without changing path length
        nor mole fractions, i.e, all quantities can be recomputed from themselves...

    Returns
    -------
    derivation: dict
        Format::

            {spectral_quantity: [list of combinations of spectral quantities needed to calculate it]}

    Examples
    --------
    to recompute a Spectrum under nonequilibrium, non optically thin case
    (note that all paths are not there yet)::

        {'abscoeff': [['absorbance']],
         'absorbance': [['transmittance_noslit'], ['abscoeff']],
         'emisscoeff': [['radiance_noslit', 'abscoeff']],
         'radiance_noslit': [['emisscoeff', 'abscoeff']],
         'transmittance_noslit': [['absorbance']]}

    a Spectrum under nonequilibrium, with optically thin case::

        {'abscoeff': [['absorbance']],
         'absorbance': [['transmittance_noslit'], ['abscoeff']],
         'emisscoeff': [['radiance_noslit']],
         'radiance_noslit': [['emisscoeff']],
         'transmittance_noslit': [['absorbance']]}

    a Spectrum under equilibrium (everything leads to everything)::

        {'abscoeff': [['absorbance'],
                      ['absorbance'],
                      ['emisscoeff'],
                      ['emissivity_noslit'],
                      ['transmittance'],
                      ['radiance'],
                      ['radiance_noslit'],
                      ['transmittance_noslit']],
         'absorbance': [['transmittance_noslit'],
                      ['abscoeff'],
                      ['abscoeff'],
                      ['emisscoeff'],
                      ['emissivity_noslit'],
                      ['transmittance'],
                      ['radiance'],
                      ['radiance_noslit'],
                      ['transmittance_noslit']],
         etc. }

    Try it ::

        import radis
        s = radis.test_spectrum()
        _build_update_graph(s)

    Outputs ::

        {'transmittance_noslit': [['absorbance'],
          ['abscoeff'],
          ['absorbance'],
          ['emisscoeff'],
          ['emissivity_noslit'],
          ['transmittance'],
          ['radiance'],
          ['radiance_noslit'],
          ['xsection']],
         'absorbance': [['transmittance_noslit'],
          ['abscoeff'],
          ['abscoeff'],
          ['emisscoeff'],
          ['emissivity_noslit'],
          ['transmittance'],
          ['radiance'],
          ['radiance_noslit'],
          ['transmittance_noslit'],
          ['xsection']],
         'abscoeff': [['absorbance'],
          ['xsection'],
          ['absorbance'],
          ['emisscoeff'],
          ['emissivity_noslit'],
          ['transmittance'],
          ['radiance'],
          ['radiance_noslit'],
          ['transmittance_noslit'],
          ['xsection']],
         'radiance_noslit': [['emisscoeff', 'abscoeff'],
          ['abscoeff'],
          ['absorbance'],
          ['emisscoeff'],
          ['emissivity_noslit'],
          ['transmittance'],
          ['radiance'],
          ['transmittance_noslit'],
          ['xsection']],
         'emisscoeff': [['radiance_noslit', 'abscoeff'],
          ['abscoeff'],
          ['absorbance'],
          ['emissivity_noslit'],
          ['transmittance'],
          ['radiance'],
          ['radiance_noslit'],
          ['transmittance_noslit'],
          ['xsection']],
         'xsection': [['abscoeff']],
         'emissivity_noslit': [['abscoeff'],
          ['absorbance'],
          ['emisscoeff'],
          ['transmittance'],
          ['radiance'],
          ['radiance_noslit'],
          ['transmittance_noslit'],
          ['xsection']]}
    """
    # Get defaults
    if path_length is None:
        path_length = "path_length" in spec.conditions
    if optically_thin is None:
        optically_thin = spec.is_optically_thin()
    if p_T_x is None:
        p_T_x = (
            "Tgas" in spec.conditions
            and "pressure_mbar" in spec.conditions
            and "mole_fraction" in spec.conditions
        )
    if equilibrium is None:
        # Read the Spectrum conditions. By default, use False.
        equilibrium = spec.conditions.get(
            "thermal_equilibrium", False
        )  # is_at_equilibrium()
    slit = (
        "slit_function" in spec.conditions
        and "slit_unit" in spec.conditions
        and "norm_by" in spec.conditions
    )

    all_usual_keys = [
        "abscoeff",
        "absorbance",
        "emisscoeff",
        "emissivity_noslit",
        "transmittance",
        "radiance",
        "radiance_noslit",
        "transmittance_noslit",
        "xsection",
    ]
    all_usual_keys_nonconvoluted = [
        k for k in all_usual_keys if k in NON_CONVOLUTED_QUANTITIES
    ]
    assert all_in(all_usual_keys, CONVOLUTED_QUANTITIES + NON_CONVOLUTED_QUANTITIES)

    # Build edges of relationships
    derivation = {  # {keys, [list of keys]}
        "transmittance_noslit": [["absorbance"]],
        "absorbance": [["transmittance_noslit"]],
    }

    def derives_from(what, *from_keys):
        """Writes that quantity ``what`` can be infered by having all
        quantities ``from_keys``

        Examples
        --------
        Radiance can be infered from emisscoeff if optically thin::

            derives_from('radiance_noslit', 'emisscoeff')

        Else abscoeff would also be needed::

            derives_from('emisscoeff', ['radiance_noslit', 'abscoeff'])
        """
        for k in from_keys:
            if isinstance(k, str):
                k = [k]
            try:
                derivation[what].append(k)
            except KeyError:
                derivation[what] = [k]

    # Build more equivalence relationships if path_length is given
    # ------------------------------------------------------------
    # TODO: complete that, list is not exhaustive yet.
    # Duplicates are removed afterwards anyway
    #
    # Note for Developers: all derives_from relationship should correspond to a
    # rescale method that was implemented. Only the developer can know that!
    # If a rescaled relationship is implemetend but not added here it wont be
    # used by the code when trying to add all quantities. If a relationship is
    # added here but not implemented it will crash during rescale (and explain why)

    # no change: case where we are rescaling without changing path length nor
    # mole fractions, i.e, all initial quantities can be recomputed from themselves...
    if no_change:
        for k in all_usual_keys:
            derives_from(k, [k])

    # Deal with case where we know path_length:
    if path_length:
        derives_from("abscoeff", ["absorbance"])
        derives_from("absorbance", ["abscoeff"])
        if optically_thin:
            derives_from("radiance_noslit", ["emisscoeff"])
            derives_from("emisscoeff", ["radiance_noslit"])
        else:
            derives_from("radiance_noslit", ["emisscoeff", "abscoeff"])
            derives_from("emisscoeff", ["radiance_noslit", "abscoeff"])

    if p_T_x:
        derives_from("abscoeff", ["xsection"])
        derives_from("xsection", ["abscoeff"])

    if slit:
        if __debug__:
            printdbg("... build_graph: slit given > convoluted keys can be recomputed")
        derives_from("radiance", ["radiance_noslit"])
        derives_from("transmittance", ["transmittance_noslit"])
        derives_from("emissivity", ["emissivity_noslit"])

    # At equilibrium, any usual spectral array can be recomputed from any other
    if equilibrium == True:
        if __debug__:
            printdbg("... build_graph: equilibrium > all keys derive from one")
        for key in all_usual_keys:
            if key in NON_CONVOLUTED_QUANTITIES and key != "xsection":
                # except from xsection, where we still need P & T
                all_but_k = [[k] for k in all_usual_keys_nonconvoluted if k != key]
                derives_from(key, *all_but_k)

    # ------------------------------------------------------------

    # if __debug__:
    #     printdbg("_build_update_graph: dependence/equivalence tree:")
    #     printdbg(derivation)

    return derivation


def get_redundant(spec):
    """Returns a dictionary of all spectral quantities in spectrum and whether
    they are redundant.

    Examples
    --------

    ::

        redundant = get_redundant(spec)

    >>> {'abscoeff': False, 'emisscoeff': True, 'absorbance': True, 'radiance_noslit': False,
         'transmittance_noslit': True}
    """

    # Get all spectral quantities that derive from existing ones
    derivation_graph = _build_update_graph(spec)

    activated = dict().fromkeys(ordered_keys, False)
    for k in spec.get_vars():
        activated[k] = True
    redundant = dict().fromkeys(ordered_keys, False)  # nothing redundant by default

    # Parse graph
    for key in ordered_keys[::-1]:  # roots
        if key in derivation_graph:
            for from_keys in derivation_graph[key]:
                if all([activated[k] and not redundant[k] for k in from_keys]):
                    redundant[key] = True  # declare redundant
                    continue
        else:
            if key not in spec.get_vars():
                del redundant[key]

    return redundant


def _path_is_complete(list_of_keys, computed_keys):
    return all([computed_keys[k] for k in list_of_keys])


def get_reachable(spec):  # , derivation_graph):
    """Get the list of all quantities that can be derived from current
    available quantities, based on given spec conditions.

    Parameters
    ----------

    spec: Spectrum
        a :class:`~radis.spectrum.spectrum.Spectrum` object

    Returns
    -------

    reachable: list
        list of quantities that can be calculated from available information

    Notes
    -----

    Algorithm::

        for all quantities, starting from the last:
            for all possible ways to compute them
                if valid, add quantity to reachable list, and restart
                else, continue
    """

    # Get inheritance links based on Spectrum conditions (equilibrium, optically thin, etc.)
    derivation_graph = _build_update_graph(spec)

    #    activated = dict().fromkeys(ordered_keys, False)
    reachable = dict().fromkeys(ordered_keys, False)
    for k in spec.get_vars():
        reachable[k] = True

    # Parse graph
    restart = True
    while restart:
        restart = False
        for key in ordered_keys[::-1]:  # roots
            if key in derivation_graph:
                # all different ways to compute this value
                for from_keys in derivation_graph[key]:
                    # if all are reachable then we can reach this new value
                    if _path_is_complete(from_keys, reachable):
                        if not reachable[key]:
                            reachable[key] = True
                            restart = True  # status changed -> restart?

    return reachable


# , derivation_graph):
def get_recompute(spec, wanted, no_change=False, true_path_length=None):
    """Get the list of all quantities that need to be recomputed to get the
    ``wanted`` quantities based on given spec conditions (does not recompute
    yet!)

    Parameters
    ----------
    spec: Spectrum
        a :class:`~radis.spectrum.spectrum.Spectrum` object
    wanted: list
        list of quantities to recompute


    Other Parameters
    ----------------
    no_change: boolean
        if True, signals that we are somehow rescaling without changing path length
        nor mole fractions, i.e, all quantities can be recomputed from themselves...
    true_path_length: boolean
        know whether the path length is given to build the equivalence graph tree.
        If None, ``path_length`` is looked up in the Spectrum condition. Default ``None``

    Returns
    -------
    recompute: list
        list of quantities needed

    Notes
    -----

    Algorithm::

        for all quantities:
            for all possible ways to compute them
                if valid, add quantity to reachable list, and restart
                else, continue
    """

    # Get inheritance links based on Spectrum conditions (equilibrium, optically thin, etc.)
    derivation_graph = _build_update_graph(
        spec, path_length=true_path_length, no_change=no_change
    )

    # Create dictionary which contains all quantities
    # that we have to recompute
    recompute = dict().fromkeys(ordered_keys, False)
    for k in spec.get_vars():  # start from all quantities we have
        recompute[k] = True
    for k in wanted:  # add all quantities we want
        recompute[k] = True
    # Add non-slit quantities if the convoluted array is wanted
    for k in CONVOLUTED_QUANTITIES:
        if k in wanted:
            recompute[k + "_noslit"] = True

    def parse_tree(recompute):
        for key in ordered_keys:
            if key in wanted:  # find a way to recompute it:
                can_recompute_key = False
                # all different ways to compute this value
                for from_keys in derivation_graph[key]:
                    # if all are recomputed already they we can recompute this new value
                    if _path_is_complete(from_keys, recompute):
                        # we can reach this quantity
                        recompute[key] = True
                        can_recompute_key = True
                        break
                if not can_recompute_key:
                    if __debug__:
                        printdbg(
                            "... get_recompute: Can't recompute {0} with current keys: {1}. Finding something else".format(
                                key, [k for k in from_keys if recompute[k]]
                            )
                        )
                    # cant recompute this quantity. Let's force recomputation
                    # of a given path. We'll arbitrary use the one will the fewer
                    # amount of not already recomputed quantities
                    score_non_recomputed_per_path = {}
                    for from_keys in derivation_graph[key]:
                        score = sum([not recompute[k] for k in from_keys])
                        # if several paths have same score the last will be chosen:
                        score_non_recomputed_per_path[score] = from_keys
                    # get the path with the minimum of values to recompute:
                    min_path = score_non_recomputed_per_path[
                        min(score_non_recomputed_per_path)
                    ]
                    # Add all these variables to the recompute list, and restart
                    for k in min_path:
                        recompute[k] = True
                        if __debug__:
                            printdbg(
                                "... get_recompute: Added new quantity to recompute list:",
                                k,
                            )
                    return recompute, True
        # reached the end with no change, no need to restart
        return recompute, False

    # Parse graph
    restart = True
    while restart:
        recompute, restart = parse_tree(recompute)

    recompute = [k for k in recompute if recompute[k]]

    if __debug__:
        printdbg("... get_recompute: List of quantities to recompute: ", recompute)

    return recompute


def update(
    spec,
    quantity="all",
    optically_thin="default",
    assume_equilibrium="default",
    verbose=True,
):
    """Calculate missing quantities that can be derived from the current
    quantities and conditions.

    e.g: if `path_length` and `emisscoeff` are given, `radiance_noslit` can be recalculated
    if `abscoeff` is also given, or without `abscoeff` in an optically thin configuration, else


    Parameters
    ----------
    spec: Spectrum
    quantity: str, ``'same'``, ``'all'``, or list of str
        name of the spectral quantity to recompute. If ``'same'``, only the quantities
        in the Spectrum are recomputed. If ``'all'``, then all quantities that can
        be derived are recomputed. Default ``'all'``. See
        :py:data:`~radis.spectrum.utils.CONVOLUTED_QUANTITIES`
        and :py:data:`~radis.spectrum.utils.NON_CONVOLUTED_QUANTITIES`
    optically_thin: True, False, or ``'default'``
        determines whether to calculate radiance with or without self-absorption.
        If ``'default'``, the value is determined from the ``'self_absorption'`` key
        in Spectrum.conditions. If not given, ``False`` is taken. Default ``'default'``
        Also updates the ``'self_absorption'`` key in conditions (creates it if
        doesnt exist)
    assume_equilibrium: boolean
        if ``True``, only absorption coefficient ``abscoeff`` is recomputed
        and all values are derived from a calculation under equilibrium,
        using Kirchoff's Law. If ``default``, use the value stored in
        Spectrum conditions[``thermal_equilibrium``], else use ``False``.
    """

    # Check inputs
    # ------------

    # Get path length
    if "path_length" in list(spec.conditions.keys()):
        path_length_cm = spec.conditions["path_length"]
        true_path_length = True
    else:
        path_length_cm = 1  # we dont know the exact path length, but some stuff can still be updated.
        true_path_length = False

    # Update thermal equilibrium
    if assume_equilibrium not in [True, False, "default"]:
        raise ValueError("assume_equilibrium must be one of True, False, 'default'")
    if assume_equilibrium == "default":
        assume_equilibrium = spec.conditions.get("thermal_equilibrium", False)
        assert assume_equilibrium in [True, False]  # prevent 'N/A' for instance

    # Update optically thin
    if optically_thin not in [True, False, "default"]:
        raise ValueError("optically_thin must be one of True, False, 'default'")
    if optically_thin == "default":
        if "self_absorption" in list(spec.conditions.keys()):
            optically_thin = not spec.conditions["self_absorption"]
        else:
            optically_thin = False
    assert optically_thin in [True, False]
    old_self_absorption = spec.conditions.get("self_absorption")
    if old_self_absorption != (not optically_thin):
        spec.conditions["self_absorption"] = not optically_thin
        if verbose:
            print(("self absorption set to:", spec.conditions["self_absorption"]))

    initial = spec.get_vars()

    # This is where everything happens:
    # ---------------------------------
    _recalculate(
        spec,
        quantity,
        # same path length (absolute value can matter)
        path_length_cm,
        path_length_cm,
        # same mole fractions (only the ratio matters)
        1,
        1,
        true_path_length=true_path_length,
        assume_equilibrium=assume_equilibrium,
        verbose=verbose,
    )

    # Output
    # ------

    # Get list of new quantities
    new_q = [k for k in spec.get_vars() if k not in initial]
    if verbose:
        if len(new_q) > 0:
            print(("{0} new quantities added: {1}".format(spec.get_name(), new_q)))

    # Final checks
    for k in new_q:
        # make sure units exist
        if not k in spec.units:
            raise ValueError("{0} added but unit is unknown".format(k))

    return spec


# Rescale functions

# ... absorption coefficient


def rescale_abscoeff(
    spec,
    rescaled,
    initial,
    old_mole_fraction,
    new_mole_fraction,
    old_path_length_cm,
    wunit,
    units,
    extra,
    true_path_length,
    assume_equilibrium,
):
    r"""

    Parameters
    ----------

    spec: Spectrum

    old_path_length: float
        path length in cm

    References
    ----------
    Scale absorption coefficient :
    .. math::

        k_2=k_1 \frac{x_2}{x_1}

    Or recompute from the absorbance:
    .. math::

        k_2=\frac{A_1}{L_1} \cdot \frac{x_2}{x_1}

    .. note::

        scaling mole fractions is not strictly valid because the lineshape changes slightly if resonant
        and non-resonant broadening coefficients are different

    Or from the transmittance :
    .. math::

        k_2=\frac{-\ln(T_1)}{L_1} \cdot \frac{x_2}{x_1}

    Or from the cross-sections :

    .. math::

        \k_2 =\sigma_1 * \frac{x p}{k_b T}  \cdot \frac{x_2}{x_1}

    """
    if __debug__:
        printdbg(
            f"recomputing `abscoeff` from initial {initial} and already rescaled {list(rescaled.keys())}, knowing `old_mole_fraction={old_mole_fraction}`, `new_mole_fraction={new_mole_fraction}`, `old_path_length_cm={old_path_length_cm}`, `true_path_length={true_path_length}`, and `assume_equilibrium={assume_equilibrium}`"
        )

    unit = None

    # case where we recomputed it already (somehow... ex: no_change signaled)
    if "abscoeff" in rescaled:
        if __debug__:
            printdbg("... rescale: abscoeff was scaled already")
        assert "abscoeff" in units
        return rescaled, units

    # First get initial abscoeff
    # ---------------------
    if "abscoeff" in initial:
        _, abscoeff_init = spec.get(
            "abscoeff", wunit=wunit, Iunit=spec.units["abscoeff"], copy=False
        )
    elif "absorbance" in initial and true_path_length:  # aka: true path_lengths given
        if __debug__:
            printdbg("... rescale: abscoeff k_1 = A_1/L_1")
        _, A = spec.get(
            "absorbance", wunit=wunit, Iunit=spec.units["absorbance"], copy=False
        )
        assert (
            units["absorbance"] == ""
        )  # TODO: allow different units, setting Iunit above?
        abscoeff_init = A / old_path_length_cm  # recalculate initial
        unit = "cm-1"
    elif "transmittance_noslit" in initial and true_path_length:
        if __debug__:
            printdbg("... rescale: abscoeff k_1 = -ln(T_1)/L_1")
        # Get abscoeff from transmittance
        _, T1 = spec.get(
            "transmittance_noslit",
            wunit=wunit,
            Iunit=spec.units["transmittance_noslit"],
            copy=False,
        )

        # We'll have a problem if the spectrum is optically thick
        b = T1 == 0  # no transmittance: optically thick mask
        if b.sum() > 0:
            msg = "Transmittance is saturated. Can't infer abscoeff. Please give absorbance"
            if "abscoeff" in extra:  # cant calculate this one but let it go
                abscoeff_init = None
                if __debug__:
                    printdbg(msg)
            else:
                raise ValueError(msg)

        # Else, let's calculate it
        abscoeff_init = -ln(T1) / old_path_length_cm  # recalculate initial
        unit = "cm-1"
    elif (
        "xsection" in initial
        and "Tgas" in spec.conditions
        and "pressure_mbar" in spec.conditions
        and "mole_fraction" in spec.conditions
    ):
        if __debug__:  # cm-1
            printdbg("... rescale: abscoeff k_2 = XS_2 * (x * p) / (k_b * T)")
        xsection = rescaled["xsection"]  # x already scaled
        pressure_Pa = spec.conditions["pressure_mbar"] * 1e2
        x = spec.conditions["mole_fraction"]
        Tgas = spec.conditions["Tgas"]  # K
        from radis.phys.constants import k_b

        abscoeff_init = xsection * (x * pressure_Pa / k_b / Tgas) * 1e-6  # cm2
        if "xsection" in spec.units:
            assert spec.units["xsection"] == "cm2"
        unit = "cm-1"

    elif "abscoeff" in extra:  # cant calculate this one but let it go
        abscoeff_init = None
    elif assume_equilibrium and "emissivity_noslit" in initial:
        raise NotImplementedError(
            "Recompute `abscoeff` from `emissivity_noslit` at equilibrium is possible but not implemented"
        )
    elif assume_equilibrium and "radiance_noslit" in initial:
        raise NotImplementedError(
            "Recompute `abscoeff` from `radiance_noslit` at equilibrium is possible but not implemented"
        )
    elif assume_equilibrium and "transmittance_noslit" in initial:
        raise NotImplementedError(
            "Recompute `abscoeff` from `transmittance_noslit` at equilibrium is possible but not implemented"
        )
    else:
        raise ValueError(
            "Can't rescale abscoeff if not all of the following are given : transmittance_noslit ({0}) ".format(
                "transmittance_noslit" in initial
            )
            + "or absorbance ({0}), ".format("absorbance" in initial)
            + "and true_path_length ({0}). ".format(true_path_length)
        )

    # Then export rescaled value
    # --------------------
    if abscoeff_init is not None:
        if __debug__:
            printdbg("... rescale: abscoeff k_2 = k_1 * (x_2/x_1)")
        abscoeff = abscoeff_init * new_mole_fraction / old_mole_fraction  # rescale
        rescaled["abscoeff"] = abscoeff
    if unit is not None:
        units["abscoeff"] = unit

    return rescaled, units


# ... all, if equilibrium and abscoeff was rescaled


def _recompute_from_abscoeff_at_equilibrium(
    spec,
    rescaled,
    wavenumber,
    Tgas,
    new_path_length_cm,
    true_path_length,
    units,
    recompute,
):
    """

    Parameters
    ----------
    rescaled: dict
        abscoeff must be rescaled

    Other Parameters
    ----------------
    recompute: list of arrays that needs to be recomputed.
    """

    def get_unit_radiance():
        """In which unit to store radiance : use the existing one by default,
        or use "mW/cm2/sr/nm if spectrum stored in wavelength; 'mW/cm2/sr/cm-1'
        if Spectrum stored in wavenumber"""
        if spec.get_waveunit() == "cm-1":
            return spec.units.get("radiance_noslit", "mW/cm2/sr/cm-1")
        else:
            return spec.units.get("radiance_noslit", "mW/cm2/sr/nm")

    def get_unit_emisscoeff(unit_radiance):
        """Basically, unit of ``unit_radiance`` divided by ``cm``"""
        if "/cm2" in unit_radiance:
            return unit_radiance.replace("/cm2", "/cm3")
        else:
            return (
                unit_radiance + "/cm"
            )  # will be simplified by Pint/Astropy.units afterwards

    assert true_path_length

    # ----------------------------------------------------------------------

    abscoeff = rescaled["abscoeff"]
    # Store:
    units["abscoeff"] = "cm-1"

    if (
        "absorbance" in recompute
        or "transmittance_noslit" in recompute
        or "emissivity_noslit" in recompute
        or "radiance_noslit" in recompute
        or "emisscoeff" in recompute
    ):
        # Calculate absorbance
        absorbance = abscoeff * new_path_length_cm
        # Store:
        rescaled["absorbance"] = absorbance
        units["absorbance"] = ""

    if (
        "transmittance_noslit" in recompute
        or "emissivity_noslit" in recompute
        or "radiance_noslit" in recompute
        or "emisscoeff" in recompute
    ):
        # Calculate transmittance
        transmittance_noslit = exp(-absorbance)
        # Store:
        rescaled["transmittance_noslit"] = transmittance_noslit
        units["transmittance_noslit"] = ""

    if (
        "emissivity_noslit" in recompute
        or "radiance_noslit" in recompute
        or "emisscoeff" in recompute
    ):
        # Calculate emissivity
        emissivity_noslit = 1 - transmittance_noslit
        # Store:
        rescaled["emissivity_noslit"] = emissivity_noslit
        units["emissivity_noslit"] = ""

    if "radiance_noslit" in recompute or "emisscoeff" in recompute:
        # Calculate radiance
        Iunit_radiance = get_unit_radiance()
        radiance_noslit = calc_radiance(
            wavenumber, emissivity_noslit, Tgas, unit=Iunit_radiance
        )
        # Store:
        rescaled["radiance_noslit"] = radiance_noslit
        units["radiance_noslit"] = Iunit_radiance

    if "emisscoeff" in recompute:
        # Calculate emission coefficient
        b = transmittance_noslit == 1  # optically thin mask
        emisscoeff = np.empty_like(abscoeff)
        emisscoeff[b] = (
            radiance_noslit[b] / new_path_length_cm
        )  # recalculate (opt thin)
        emisscoeff[~b] = (
            radiance_noslit[~b] / (1 - transmittance_noslit[~b]) * abscoeff[~b]
        )  # recalculate (non opt thin)
        # Store:
        rescaled["emisscoeff"] = emisscoeff
        units["emisscoeff"] = get_unit_emisscoeff(units["radiance_noslit"])

    return rescaled, units


# ... emission coefficient


def rescale_emisscoeff(
    spec,
    rescaled,
    initial,
    old_mole_fraction,
    new_mole_fraction,
    old_path_length_cm,
    optically_thin,
    wunit,
    units,
    extra,
    true_path_length,
):
    r"""

    Parameters
    ----------
    spec: Spectrum


    References
    ----------
    If optically thin, compute from the radiance :

    .. math::

        j_2=\frac{I_1}{L_1} \cdot \frac{x_2}{x_1}


    In the general case, compute from the radiance and the absorption coefficient :

    .. math::

        j_2=\frac{k_1 I_1}{1-\operatorname{exp}\left(-k_1 L_1\right)} \cdot \frac{x_2}{x_1}

    .. note::

        scaling path length is always valid ; scaling mole fractions is not
        strictly valid because the lineshape changes slightly if resonant
        and non-resonant broadening coefficients are different

    or from the radiance and the transmittance:

    .. math::

        j_2=\frac{k_1 I_1}{1-T_1} \frac{x_2}{x_1}


    """

    unit = None

    def get_emisscoeff_unit(unit_radiance):
        """Basically, units of ``unit_radiance`` divided by ``cm``"""
        if "/cm2" in unit_radiance:
            return unit_radiance.replace("/cm2", "/cm3")
        else:
            return (
                unit_radiance + "/cm"
            )  # will be simplified by Pint/Astropy.units afterwards

    # case where we recomputed it already (somehow... ex: no_change signaled)
    if "emisscoeff" in rescaled:
        if __debug__:
            printdbg("... rescale: emisscoeff was scaled already")
        assert "emisscoeff" in units
        return rescaled, units

    # Firt get initial emisscoeff j1
    # -------------------

    if "emisscoeff" in initial:
        if __debug__:
            printdbg("... rescale: emisscoeff j1 = j1")
        _, emisscoeff_init = spec.get(
            "emisscoeff", wunit=wunit, Iunit=units["emisscoeff"], copy=False
        )

    elif "radiance_noslit" in initial and true_path_length and optically_thin:
        if __debug__:
            printdbg("... rescale: emisscoeff j_1 = I_1/L_1")
        _, I = spec.get(
            "radiance_noslit", wunit=wunit, Iunit=units["radiance_noslit"], copy=False
        )
        emisscoeff_init = I / old_path_length_cm  # recalculate initial
        unit = get_emisscoeff_unit(units["radiance_noslit"])

    elif "radiance_noslit" in initial and true_path_length and "abscoeff" in initial:
        if __debug__:
            printdbg("... rescale: emisscoeff j_1 = k_1*I_1/(1-exp(-k_1*L_1))")
        # get emisscoeff from (initial) abscoeff and (initial) radiance
        _, I = spec.get("radiance_noslit", wunit=wunit, Iunit=units["radiance_noslit"])
        _, k = spec.get("abscoeff", wunit=wunit, Iunit=units["abscoeff"], copy=False)

        # Recalculate in the optically thin range (T=1) and elsewhere
        b = k == 0  # optically thin mask
        emisscoeff_init = np.empty_like(k)
        # ... optically thin case
        # recalculate (opt thin)
        emisscoeff_init[b] = I[b] / old_path_length_cm

        # ... non optically thin case:
        # ... recalculate transmittance from abscoeff
        T_b = exp(-k[~b] * old_path_length_cm)  # recalculate initial
        # ... and solve the RTE on an homogeneous slab
        # recalculate (non opt thin)
        emisscoeff_init[~b] = k[~b] * I[~b] / (1 - T_b)
        unit = get_emisscoeff_unit(unit_radiance=units["radiance_noslit"])

    elif (
        "radiance_noslit" in initial
        and true_path_length
        and "transmittance_noslit" in initial
    ):
        if __debug__:
            printdbg("... rescale: emisscoeff j_1 = k_1*I_1/(1-T_1)")
        # get emisscoeff from (initial) transmittance and (initial) radiance
        _, I = spec.get("radiance_noslit", wunit=wunit, Iunit=units["radiance_noslit"])
        _, T = spec.get(
            "transmittance_noslit", wunit=wunit, Iunit=units["transmittance_noslit"]
        )

        # Recalculate in the optically thin range (T=1) and elsewhere
        b = T == 1  # optically thin mask
        emisscoeff_init = np.empty_like(T)
        # ... optically thin case
        # recalculate (opt thin)
        emisscoeff_init[b] = I[b] / old_path_length_cm

        # ... non optically thin case:
        # ... recalculate abscoeff from transmittance
        T_b = T[~b]
        k_b = -ln(T_b) / old_path_length_cm  # cm-1
        # ... and solve the RTE on an homogeneous slab
        # recalculate (non opt thin)
        emisscoeff_init[~b] = k_b * I[~b] / (1 - T_b)
        unit = get_emisscoeff_unit(unit_radiance=units["radiance_noslit"])

    else:
        if optically_thin:
            msg = "Can't calculate emisscoeff if path_length ({0})".format(
                true_path_length
            ) + "and initial radiance_noslit ({0}) are not all given".format(
                "radiance_noslit" in initial
            )
            if "emisscoeff" in extra:  # cant calculate this one but let it go
                emisscoeff_init = None
                if __debug__:
                    printdbg(msg)
            else:
                raise ValueError(msg)
        else:
            msg = (
                "Trying to get the emission coefficient (emisscoeff) for a non-optically "
                + "thin column. Among the following required quantities, not all of them are given : path_length ({0}), radiance_noslit ({1}) ".format(
                    true_path_length, "radiance_noslit" in initial
                )
                + "and abscoeff ({0}). ".format("abscoeff" in initial)
                + "See known Spectrum conditions with "
                + "print(Spectrum)"
            )
            if "emisscoeff" in extra:  # cant calculate this one but let it go
                emisscoeff_init = None
                if __debug__:
                    printdbg(msg)
            else:
                raise ValueError(msg)

    # Then rescale and export
    # -----------------------
    if emisscoeff_init is not None:
        if __debug__:
            printdbg("... rescale: emisscoeff j_2 = j_1 * (x_2/x_1)")
        # Now rescale for mole fractions
        emisscoeff = (
            emisscoeff_init * new_mole_fraction / old_mole_fraction
        )  # rescale x
        rescaled["emisscoeff"] = emisscoeff
    if unit is not None:
        units["emisscoeff"] = unit

    return rescaled, units


# ... absorbance


def rescale_absorbance(
    spec,
    rescaled,
    initial,
    old_mole_fraction,
    new_mole_fraction,
    old_path_length_cm,
    new_path_length_cm,
    waveunit,
    rescaled_units,
    extra,
    true_path_length,
):
    """

    Parameters
    ----------
    spec: Spectrum

    References
    ----------
    Rescale the absorbance:

    .. math::

        A_2=A_1 \frac{x_2}{x_1} \frac{L_2}{L_1}

    .. note::

        scaling path length is always valid ; scaling mole fractions is not
        strictly valid because the lineshape changes slightly if resonant
        and non-resonant broadening coefficients are different

    Or compute from the absorption coefficient :

    .. math::

        A_2=k_2 L_2

    Or from the transmittance :

    .. math::

        A_2=-\ln(T1) \frac{x_2}{x_1} \frac{L_2}{L_1}

    See Also
    --------
    :py:attr:`~radis.spectrum.rescale.rescale_abscoeff`,
    :py:attr:`~radis.spectrum.rescale.rescale_transmittance_noslit`,

    """

    unit = None

    # case where we recomputed it already (somehow... ex: no_change signaled)
    if "absorbance" in rescaled:
        if __debug__:
            printdbg("... rescale: absorbance was scaled already")
        assert "absorbance" in rescaled_units
        return rescaled, rescaled_units

    # Get scaled absorbance directly A1
    # ---------------------------------

    if "absorbance" in initial:
        if __debug__:
            printdbg("... rescale: absorbance A_2 = A_1*(x_2/x_1)*(L_2/L_1)")
        _, absorbance = spec.get(
            "absorbance", wunit=waveunit, Iunit=spec.units["absorbance"], copy=True
        )
        assert spec.units["absorbance"] == ""
        absorbance *= new_mole_fraction / old_mole_fraction  # rescale x
        absorbance *= new_path_length_cm / old_path_length_cm  # rescale L
        unit = ""
    elif "abscoeff" in rescaled and true_path_length:  # in cm
        if __debug__:
            printdbg("... rescale: absorbance A_2 = k_2*L_2")
        abscoeff = rescaled["abscoeff"]  # x already scaled
        absorbance = abscoeff * new_path_length_cm  # calculate L
        assert rescaled_units["abscoeff"] == "cm-1"
        unit = ""
    elif "transmittance_noslit" in initial and true_path_length:
        if __debug__:
            printdbg("... rescale: absorbance A_2 = -ln(T1)*(x_2/x_1)*(L_2/L_1)")
        # Get absorbance from transmittance
        _, T1 = spec.get(
            "transmittance_noslit",
            wunit=waveunit,
            Iunit=spec.units["transmittance_noslit"],
            copy=False,
        )

        # We'll have a problem if the spectrum is optically thick
        b = T1 == 0  # no transmittance: optically thick mask
        if b.sum() > 0:
            msg = "Transmittance is saturated. Can't calculate absorbance"
            if "absorbance" in extra:  # cant calculate this one but let it go
                absorbance = None
                if __debug__:
                    printdbg(msg)
            else:
                raise ValueError(msg)

        # Else, let's calculate it
        absorbance = -ln(T1)
        absorbance *= new_mole_fraction / old_mole_fraction  # rescale x
        absorbance *= new_path_length_cm / old_path_length_cm  # rescale L
        unit = ""
    else:
        msg = (
            "Cant recalculate absorbance if scaled absoeff "
            + "({0}) and true path_length ({1}) are not given".format(
                "abscoeff" in rescaled, true_path_length
            )
        )
        if "absorbance" in extra:  # cant calculate this one but let it go
            absorbance = None
            if __debug__:
                printdbg(msg)
        else:
            raise ValueError(msg)

    # Then rescale and export
    # -----------------------

    # Export rescaled value
    if absorbance is not None:
        rescaled["absorbance"] = absorbance
    if rescaled_units is not None:
        rescaled_units["absorbance"] = unit

    return rescaled, rescaled_units


# ... transmittance
def rescale_transmittance_noslit(
    spec,
    rescaled,
    initial,
    old_mole_fraction,
    new_mole_fraction,
    old_path_length_cm,
    new_path_length_cm,
    waveunit,
    units,
    extra,
    true_path_length,
):
    r"""

    Parameters
    ----------
    spec: Spectrum


    References
    ----------
    Compute from the absorbance :

    .. math::

        \tau_2=\operatorname{exp}\left(-A_2\right)

    Or from the absorption coefficient and the new path length :

    .. math::

        \tau_2=\operatorname{exp}\left(-k_2 L_2\right)

    Or from the previous transmittance :

    .. math::

        \tau_2=\operatorname{exp}\left(\ln(T_1) \frac{x_2}{x_1} \frac{L_2}{L_1}\right)

    .. note::

        scaling path length is always valid ; scaling mole fractions is not
        strictly valid because the lineshape changes slightly if resonant
        and non-resonant broadening coefficients are different

    See Also
    --------
    :py:attr:`~radis.spectrum.rescale.rescale_abscoeff`,
    :py:attr:`~radis.spectrum.rescale.rescale_absorbance`,
    :py:attr:`~radis.spectrum.rescale.rescale_transmittance_noslit`,
    """

    unit = None

    def get_transmittance_unit():
        return ""

    # case where we recomputed it already (somehow... ex: no_change signaled)
    if "transmittance_noslit" in rescaled:
        if __debug__:
            printdbg("... rescale: transmittance_noslit was scaled already")
        assert "transmittance_noslit" in units
        return rescaled, units

    # Calculate rescaled value directly
    # ---------------------------------

    # Rescale
    if "absorbance" in rescaled:
        if __debug__:
            printdbg("... rescale: transmittance_noslit T_2 = exp(-A_2)")
        absorbance = rescaled["absorbance"]  # x and L already scaled
        transmittance_noslit = exp(-absorbance)  # recalculate
        unit = get_transmittance_unit()
    elif "abscoeff" in rescaled and true_path_length:
        if __debug__:
            printdbg("... rescale: transmittance_noslit T_2 = exp(-k_2*L_2)")
        abscoeff = rescaled["abscoeff"]  # x already scaled
        assert units["abscoeff"] == "cm-1"
        absorbance = abscoeff * new_path_length_cm  # calculate
        transmittance_noslit = exp(-absorbance)  # recalculate
        unit = get_transmittance_unit()
    elif "transmittance_noslit" in initial:
        if __debug__:
            printdbg(
                "... rescale: transmittance_noslit T_2 = "
                + "exp( ln(T_1) * (x_2/x_1) * (L_2/L_1))"
            )
        # get transmittance from initial transmittance
        _, T1 = spec.get(
            "transmittance_noslit",
            wunit=waveunit,
            Iunit=units["transmittance_noslit"],
            copy=False,
        )

        # We'll have a problem if the spectrum is optically thick
        b = T1 == 0  # optically thick mask
        if b.sum() > 0 and (
            new_mole_fraction < old_mole_fraction
            or new_path_length_cm < old_path_length_cm
        ):
            # decreasing mole fractions/ path length could increase the transmittance
            # but this information was lost in the saturation
            msg = "Transmittance is saturated. Cant rescale. Please give absorbance"
            if "transmittance_noslit" in extra:  # cant calculate this one but let it go
                transmittance_noslit = None
                if __debug__:
                    printdbg(msg)
            else:
                raise ValueError(msg)
        # Else, just get absorbance
        absorbance = -ln(T1)
        absorbance *= new_mole_fraction / old_mole_fraction  # rescale x
        absorbance *= new_path_length_cm / old_path_length_cm  # rescale L
        transmittance_noslit = exp(-absorbance)
        unit = get_transmittance_unit()
    else:
        msg = "Missing data to rescale transmittance. Expected scaled absorbance ({0})".format(
            "absorbance" in rescaled
        )
        #        +' or scaled abscoeff ({0}) and true_path_length ({1})'.format(
        #                'abscoeff' in rescaled, true_path_length)
        if "transmittance_noslit" in extra:  # cant calculate this one but let it go
            transmittance_noslit = None
            if __debug__:
                printdbg(msg)
        else:
            raise ValueError(msg)

    # Export rescaled value
    if transmittance_noslit is not None:
        rescaled["transmittance_noslit"] = transmittance_noslit
    if unit is not None:
        units["transmittance_noslit"] = unit

    return rescaled, units


def rescale_transmittance(
    spec,
    rescaled,
    initial,
    old_mole_fraction,
    new_mole_fraction,
    old_path_length_cm,
    new_path_length_cm,
    waveunit,
    units,
    extra,
):
    """

    Parameters
    ----------
    spec: Spectrum
    """

    #    unit = None
    apply_slit = False
    #    def get_unit():
    #        return '1'

    # case where we recomputed it already (somehow... ex: no_change signaled)
    if "transmittance" in rescaled:
        if __debug__:
            printdbg("... rescale: transmittance was scaled already")
        assert "transmittance" in units
        return rescaled, units, apply_slit

    if "transmittance_noslit" in rescaled:
        apply_slit = True
    else:
        raise NotImplementedError("rescale transmittance not implemented")

    return rescaled, units, apply_slit


# ... radiance_noslit


def rescale_radiance_noslit(
    spec,
    rescaled,
    initial,
    old_mole_fraction,
    new_mole_fraction,
    old_path_length_cm,
    new_path_length_cm,
    optically_thin,
    waveunit,
    rescaled_units,
    extra,
    true_path_length,
):
    r"""

    Parameters
    ----------
    spec: Spectrum


    References
    ----------
    If optically thin, calculate from the emission coefficient :

    .. math::

        I_2=j_2 L_2

    or scale the existing radiance :

    .. math::

        I_2=\frac{\frac{I_1 x_2}{x_1} L_2}{L_1}

    .. note::

        scaling path length is always valid ; scaling mole fractions is not
        strictly valid because the lineshape changes slightly if resonant
        and non-resonant broadening coefficients are different

    If not optically thin, recompute from the transmittance and absorption coefficient
    (analytical output of the 1D-homogenous column emission without scattering) :

    .. math::

        I_2=\frac{j_2 \left(1-T_2\right)}{k_2}

    or from the emission coefficient and absorption coefficient:

    .. math::

        I_2=\frac{j_2 \left(1-\operatorname{exp}\left(-k_2 L_2\right)\right)}{k_2}


    See Also
    --------
    :py:attr:`~radis.spectrum.rescale.rescale_emisscoeff`,
    :py:attr:`~radis.spectrum.rescale.rescale_abscoeff`,
    :py:attr:`~radis.spectrum.rescale.rescale_radiance_noslit`,

    """
    if __debug__:
        printdbg(
            f"recomputing `radiance_noslit` from initial {initial} and already rescaled {list(rescaled.keys())}, knowing `old_mole_fraction={old_mole_fraction}`, `new_mole_fraction={new_mole_fraction}`, `old_path_length_cm={old_path_length_cm}`, `new_path_length_cm={new_path_length_cm}`, `true_path_length={true_path_length}`"
        )

    unit = None

    def get_radiance_unit(unit_emisscoeff):
        """get radiance_noslit unit from emisscoeff unit.
        Basically just multiply by ``cm``"""
        if "/cm3" in unit_emisscoeff:
            return unit_emisscoeff.replace("/cm3", "/cm2")
        else:
            return unit_emisscoeff + "*cm"

    # case where we recomputed it already (somehow... ex: no_change signaled)
    if "radiance_noslit" in rescaled:
        if __debug__:
            printdbg("... rescale: radiance_noslit was scaled already")
        assert "radiance_noslit" in rescaled_units
        return rescaled, rescaled_units

    # Rescale!
    if "emisscoeff" in rescaled and true_path_length and optically_thin:
        if __debug__:
            printdbg(
                "... rescale: radiance_noslit I_2 = j_2 * L_2 " + "(optically thin)"
            )
        emisscoeff = rescaled["emisscoeff"]  # x already scaled
        radiance_noslit = emisscoeff * new_path_length_cm  # recalculate L
        unit = get_radiance_unit(rescaled_units["emisscoeff"])

    elif (
        "emisscoeff" in rescaled
        and "transmittance_noslit" in rescaled
        and "abscoeff" in rescaled
        and true_path_length
        and not optically_thin
    ):  # not optically thin
        if __debug__:
            printdbg("... rescale: radiance_noslit I_2 = j_2*(1-T_2)/k_2")
        emisscoeff = rescaled["emisscoeff"]  # x already scaled
        abscoeff = rescaled["abscoeff"]  # x already scaled
        assert rescaled_units["abscoeff"] == "cm-1"
        # mole_fraction, path_length already scaled
        transmittance_noslit = rescaled["transmittance_noslit"]
        b = transmittance_noslit == 1  # optically thin mask
        radiance_noslit = np.empty_like(emisscoeff)  # calculate L
        radiance_noslit[~b] = (
            emisscoeff[~b] / abscoeff[~b] * (1 - transmittance_noslit[~b])
        )
        radiance_noslit[b] = emisscoeff[b] * new_path_length_cm  # optically thin limit
        unit = get_radiance_unit(unit_emisscoeff=rescaled_units["emisscoeff"])

    elif (
        "emisscoeff" in rescaled
        and "abscoeff" in rescaled
        and true_path_length
        and not optically_thin
    ):  # not optically thin
        if __debug__:
            printdbg("... rescale: radiance_noslit I_2 = j_2*(1-exp(-k_2*L_2))/k_2")
        emisscoeff = rescaled["emisscoeff"]  # x already scaled
        abscoeff = rescaled["abscoeff"]  # x already scaled
        assert rescaled_units["abscoeff"] == "cm-1"
        b = abscoeff == 0  # optically thin mask
        radiance_noslit = np.empty_like(emisscoeff)  # calculate
        radiance_noslit[~b] = (
            emisscoeff[~b]
            / abscoeff[~b]
            * (1 - exp(-abscoeff[~b] * new_path_length_cm))
        )
        radiance_noslit[b] = emisscoeff[b] * new_path_length_cm  # optically thin limit
        unit = get_radiance_unit(unit_emisscoeff=rescaled_units["emisscoeff"])

    elif "radiance_noslit" in initial and optically_thin:
        if __debug__:
            printdbg(
                "... rescale: radiance_noslit I_2 = I_1*x_2/x_1*L_2/L_1 "
                + "(optically thin)"
            )
        _, radiance_noslit = spec.get(
            "radiance_noslit",
            wunit=waveunit,
            Iunit=rescaled_units["radiance_noslit"],
            copy=True,
        )
        radiance_noslit *= new_mole_fraction / old_mole_fraction  # rescale
        radiance_noslit *= new_path_length_cm / old_path_length_cm  # rescale
        unit = rescaled_units["radiance_noslit"]

    else:
        if optically_thin:
            msg = (
                "Missing data to rescale radiance_noslit in "
                + "optically thin mode. You need at least initial "
                + "radiance_noslit ({0}), or scaled emission coefficient ({1}) ".format(
                    "radiance_noslit" in initial, "emisscoeff" in rescaled
                )
                + "and true path length ({0}).".format(true_path_length)
            )
            if "radiance_noslit" in extra:  # cant calculate this one but let it go
                radiance_noslit = None
                if __debug__:
                    printdbg(msg)
            else:
                raise ValueError(msg)
        else:
            msg = (
                "Missing data to recalculate radiance_noslit for a non-optically thin column. You need at least "
                + "scaled emisscoeff ({0}), scaled transmittance_noslit ({1}), ".format(
                    "emisscoeff" in rescaled, "transmittance_noslit" in rescaled
                )
                + "scaled abscoeff ({0}) and true_path_length ({1}). ".format(
                    "abscoeff" in rescaled, true_path_length
                )
            )
            if "radiance_noslit" in extra:
                radiance_noslit = None
                if __debug__:
                    printdbg(msg)
            else:
                raise ValueError(msg)

    # Export rescaled value
    if radiance_noslit is not None:
        rescaled["radiance_noslit"] = radiance_noslit
    if unit is not None:
        rescaled_units["radiance_noslit"] = unit

    return rescaled, rescaled_units


# ... radiance


def rescale_radiance(
    spec,
    rescaled,
    initial,
    old_mole_fraction,
    new_mole_fraction,
    old_path_length_cm,
    new_path_length_cm,
    optically_thin,
    waveunit,
    rescaled_units,
    extra,
    true_path_length,
):
    """

    Parameters
    ----------

    spec: Spectrum
    """

    apply_slit = False
    #    unit = None
    #    def get_unit(unit_emisscoeff):
    #        if '/cm3' in unit_emisscoeff:
    #            return unit_emisscoeff.replace('/cm3', '/cm2')
    #        else:
    #            return unit_emisscoeff + '*cm'

    # case where we recomputed it already (somehow... ex: no_change signaled)
    if "radiance" in rescaled:
        if __debug__:
            printdbg("... rescale: radiance was scaled already")
        assert "radiance" in rescaled_units
        return rescaled, rescaled_units, apply_slit

    # Rescale!
    if "radiance_noslit" in rescaled:
        apply_slit = True
    else:
        raise NotImplementedError("rescale radiance not implemented yet")

    return rescaled, rescaled_units, apply_slit


# ... emissivity_noslit


def rescale_emissivity_noslit(spec, rescaled, rescaled_units, extra, true_path_length):
    r"""

    Parameters
    ----------

    spec: Spectrum

    References
    ----------

    Compute from the transmittance :

    .. math::

        \epsilon_2 = 1 - \tau_2

    See Also
    --------
    :py:attr:`~radis.spectrum.rescale.rescale_transmittance_noslit`,
    """

    # case where we recomputed it already (somehow... ex: no_change signaled)
    # -------------------------
    if "emissivity_noslit" in rescaled:
        if __debug__:
            printdbg("... rescale: emissivity_noslit was scaled already")
        assert "emissivity_noslit" in rescaled_units
        return rescaled, rescaled_units

    # Or: rescale!
    # -------------
    if "transmittance_noslit" in rescaled:
        if __debug__:
            printdbg("... rescale: emissivity_noslit e_2 = 1 - T_2")
        # transmittivity already scaled
        T2 = rescaled["transmittance_noslit"]
        emissivity_noslit = 1 - T2  # recalculate
    else:
        msg = "transmittance_noslit needed to recompute emissivity_noslit"
        if "emissivity_noslit" in extra:  # cant calculate this one but let it go
            emissivity_noslit = None
            if __debug__:
                printdbg(msg)
        else:
            raise ValueError(msg)

    # Export rescaled value
    if emissivity_noslit is not None:
        rescaled["emissivity_noslit"] = emissivity_noslit
        rescaled_units["emissivity_noslit"] = ""

    return rescaled, rescaled_units


# ... cross sections


def rescale_xsection(
    spec,
    rescaled,
    initial,
    old_mole_fraction,
    new_mole_fraction,
    old_path_length_cm,
    new_path_length_cm,
    waveunit,
    rescaled_units,
    extra,
    true_path_length,
):
    r"""Rescale cross-sections

    Parameters
    ----------
    spec: Spectrum
    old_path_length: float
        path length in cm

    References
    ----------
    Either scale the cross-section directly :

    .. math::

        \sigma_2=\sigma_1

    (by definition, cross-sections are unchanged)

    or recompute from scaled absorption coefficient:

    .. math::

        \sigma_2=k_2 \frac{k_b T}{x p}

    With ``p`` the total pressure and ``x`` the species mole fraction.


    See Also
    --------
    :py:attr:`~radis.spectrum.rescale.rescale_abscoeff`,
    """

    unit = None

    # case where we recomputed it already (somehow... ex: no_change signaled)
    if "xsection" in rescaled:
        if __debug__:
            printdbg("... rescale: xsection was scaled already")
        assert "xsection" in rescaled_units
        return rescaled, rescaled_units

    # Get scaled xsection directly from XS1
    # -------------------------------------

    if "xsection" in initial:
        # cross-sections are unchanged
        if __debug__:
            printdbg("... rescale: xsection XS_2 = XS_1")
        _, xsection = spec.get(
            "xsection", wunit=waveunit, Iunit=rescaled_units["xsection"]
        )
        unit = rescaled_units["xsection"]  # units unchanged
    elif (
        "abscoeff" in rescaled
        and "Tgas" in spec.conditions
        and "pressure_mbar" in spec.conditions
        and "mole_fraction" in spec.conditions
    ):
        if __debug__:  # cm-1
            printdbg("... rescale: xsection XS_2 = k_2 * (k_b * T / x / p)")
        abscoeff = rescaled["abscoeff"]  # x already scaled
        assert spec.units["abscoeff"] == "cm-1"

        xsection = abscoeff2xsection(
            abscoeff_cm1=abscoeff,
            Tgas_K=spec.conditions["Tgas"],
            mole_fraction=spec.conditions["mole_fraction"],
            pressure_Pa=spec.conditions["pressure_mbar"] * 1e2,
        )
        unit = "cm2"
    else:
        msg = (
            "Cant recalculate xsection if not all these quantities are given : "
            + "scaled abscoeff ({0}), Tgas ({1}), pressure_mbar ({2})".format(
                "abscoeff" in rescaled,
                "Tgas" in spec.conditions,
                "pressure_mbar" in spec.conditions,
            )
        )
        if "xsection" in extra:  # cant calculate this one but let it go
            xsection = None
            if __debug__:
                printdbg(msg)
        else:
            raise ValueError(msg)

    # Then rescale and export
    # -----------------------

    # Export rescaled value
    if xsection is not None:
        rescaled["xsection"] = xsection
    if unit is not None:
        rescaled_units["xsection"] = unit

    return rescaled, rescaled_units


def _recalculate(
    spec,
    quantity,
    new_path_length_cm,
    old_path_length_cm,
    new_mole_fraction,
    old_mole_fraction,
    true_path_length=True,
    verbose=True,
    assume_equilibrium=False,
):
    """General function to recalculate missing quantities. Used in
    :func:`~radis.spectrum.rescale.rescale_path_length`, :func:`~radis.spectrum.rescale.rescale_mole_fraction`
    and :func:`~radis.spectrum.rescale.update`.

    Determines with spectral quantities should be recomputed, then scales
    them solving the Radiative Transfer Equation on a 1D-homogeneous column
    without scattering.


    Parameters
    ----------
    spec: Spectrum
        the Spectrum object to recompute
    quantity: str, ``'same'``, ``'all'``, or list of str
        name of the spectral quantity to recompute. If ``'same'``, only the quantities
        in the Spectrum are recomputed. If ``'all'``, then all quantities that can
        be derived are recomputed. See :py:data:`~radis.spectrum.utils.CONVOLUTED_QUANTITIES`
        and :py:data:`~radis.spectrum.utils.NON_CONVOLUTED_QUANTITIES`
    true_path_length: boolean
        if ``False``, only relative rescaling (new/old) is allowed. For instance,
        when you dont know the true path_lenth, rescaling absorbance
        with *= new_length/old_length is fine, but abscoeff*new_length is not
        Default ``True``

    Other Parameters
    ----------------
    assume_equilibrium: boolean
        if ``True``, only absorption coefficient ``abscoeff`` is recomputed
        and all values are derived from a calculation under equilibrium,
        using Kirchoff's Law. Default ``False``
    """

    optically_thin = spec.is_optically_thin()
    initial = spec.get_vars()  # quantities initially in spectrum
    if __debug__:
        printdbg("... rescale: optically_thin: {0}".format(optically_thin))
        printdbg("... rescale: initial quantities: {0}".format(initial))

    # Check that inputs are valid names
    _check_quantity = quantity
    if isinstance(quantity, list):
        _check_quantity = quantity
    else:
        _check_quantity = [quantity]
    for k in _check_quantity:
        try:
            assert k in CONVOLUTED_QUANTITIES + NON_CONVOLUTED_QUANTITIES + [
                "all",
                "same",
            ]
        except AssertionError:
            raise ValueError(
                f"Unexpected spectral array `{k}`. Expected one of {CONVOLUTED_QUANTITIES+NON_CONVOLUTED_QUANTITIES}"
            )
    # ... make sure we're not trying to rescale a Spectrum that has non scalable
    # ... quantities
    if any_in(initial, non_rescalable_keys):
        raise NotImplementedError(
            "Trying to rescale a Spectrum that has non scalable "
            + "quantities: {0}".format([k for k in initial if k in non_rescalable_keys])
            + "Remove them manually, or implement the scaling method."
        )

    # Choose which values to recompute (and store them in the list `wanted`)
    # ----------
    if quantity == "all":  # quantities to recompute
        wanted = list(initial)  # start from the one we have (that also makes)
        # sure we dont delete anyone with the final validity
        # check "everyone is here"
        greedy = True
    elif quantity == "same":
        wanted = list(initial)
        greedy = False
    elif isinstance(quantity, str):
        wanted = [quantity]
        greedy = False
    elif isinstance(quantity, list):
        wanted = quantity
        greedy = False
    else:
        raise ValueError(
            "unexpected type for quantity: expected str, got "
            + "{0} ({1})".format(quantity, type(quantity))
        )
    rescaled = {}  # quantities rescaled

    # in greedy mode ('all'), choose to recompute all parameters that we can
    extra = []
    if greedy:
        # ... let's be greedy: recompute all possible quantities. The list of
        # all spectral quantities is calculated by parsing a tree in get_reachable
        reachable = get_reachable(spec)
        extra = [k for k, v in reachable.items() if v]
    wanted = set(wanted + extra)

    # There are two cases: either we are actually rescaling to another length /
    # mole fraction, or we are just updating() without changing length / mole fraction
    no_change = (
        new_mole_fraction == old_mole_fraction
        and new_path_length_cm == old_path_length_cm
    )

    # Quickly stop if no change
    if no_change and all_in(wanted, initial):
        if __debug__:
            printdbg("... rescale: no change")
        # Stop here
        return

    # list of quantities that are needed to recompute what we want
    # ... (we're just analysing how to compute them here, the actual calculation
    # ... will be done later)
    try:
        recompute = get_recompute(
            spec, wanted, no_change, true_path_length=true_path_length
        )
    except KeyError as err:
        import sys

        print(sys.exc_info())
        raise KeyError(
            "Error in get_recompute (see above). Quantity `{0}` cannot be recomputed ".format(
                err.args[0]
            )
            + "from available quantities in Spectrum ({0}) with ".format(
                spec.get_vars()
            )
            + " conditions: optically thin ({0}), true_path_length ({1}), thermal equilibrium ({2})".format(
                optically_thin, true_path_length, assume_equilibrium
            )
            + ". Check how your equivalence tree is built: see rescale._build_update_graph()"
        )
    if assume_equilibrium == True:
        recompute.append("abscoeff")
    recompute = set(recompute)  # remove duplicates

    # Get units
    rescaled_units = spec.units.copy()

    # Recompute!
    # ----------
    waveunit = spec.get_waveunit()  # keep all quantities in same waveunit
    apply_slit = False  # if True at the end re-apply slit

    # If no_change, just set everyone as rescaled already
    if no_change:
        for k in initial:
            try:
                spec.units[k]
            except KeyError:
                raise KeyError(f"Units for {k} must be defined")
            rescaled[k] = spec.get(k, wunit=spec.get_waveunit(), Iunit=spec.units[k])[
                1
            ]  # note: creates a copy

    # Start with abscoeff

    if "abscoeff" in recompute:
        rescaled, rescaled_units = rescale_abscoeff(
            spec,
            rescaled,
            initial,  # Todo: remove rescaled = ... Dict is mutable no?
            old_mole_fraction,
            new_mole_fraction,
            old_path_length_cm,
            waveunit,
            rescaled_units,
            extra,
            true_path_length,
            assume_equilibrium,
        )

    if (
        assume_equilibrium
        and "Tgas" in spec.conditions
        and spec.conditions["Tgas"] != "N/A"
    ):
        # ... (dont forget Python will stop at evaluating the 2nd expression if False)
        assert "abscoeff" in rescaled
        if not spec.is_at_equilibrium():
            warn(
                "Rescaling with equilibrium assumption but Spectrum {0} does ".format(
                    spec.get_name()
                    + "not look at equilibrium. Check spectrum conditions"
                )
            )
        wavenumber = spec.get_wavenumber()
        Tgas = spec.conditions["Tgas"]
        rescaled, rescaled_units = _recompute_from_abscoeff_at_equilibrium(
            spec,
            rescaled,
            wavenumber,
            Tgas,
            new_path_length_cm,
            true_path_length,
            rescaled_units,
            recompute,
        )
        apply_slit = "radiance" in recompute or "transmittance" in recompute

    else:

        if "emisscoeff" in recompute:
            rescaled, rescaled_units = rescale_emisscoeff(
                spec,
                rescaled,
                initial,
                old_mole_fraction,
                new_mole_fraction,
                old_path_length_cm,
                optically_thin,
                waveunit,
                rescaled_units,
                extra,
                true_path_length,
            )

        if "absorbance" in recompute:
            rescaled, rescaled_units = rescale_absorbance(
                spec,
                rescaled,
                initial,
                old_mole_fraction,
                new_mole_fraction,
                old_path_length_cm,
                new_path_length_cm,
                waveunit,
                rescaled_units,
                extra,
                true_path_length,
            )

        if "transmittance_noslit" in recompute:
            rescaled, rescaled_units = rescale_transmittance_noslit(
                spec,
                rescaled,
                initial,
                old_mole_fraction,
                new_mole_fraction,
                old_path_length_cm,
                new_path_length_cm,
                waveunit,
                rescaled_units,
                extra,
                true_path_length,
            )

        if "radiance_noslit" in recompute:
            rescaled, rescaled_units = rescale_radiance_noslit(
                spec,
                rescaled,
                initial,
                old_mole_fraction,
                new_mole_fraction,
                old_path_length_cm,
                new_path_length_cm,
                optically_thin,
                waveunit,
                rescaled_units,
                extra,
                true_path_length,
            )

        if "emissivity_noslit" in recompute:
            rescaled, rescaled_units = rescale_emissivity_noslit(
                spec, rescaled, rescaled_units, extra, true_path_length
            )

        if "radiance" in recompute:
            rescaled, rescaled_units, slit_needed = rescale_radiance(
                spec,
                rescaled,
                initial,
                old_mole_fraction,
                new_mole_fraction,
                old_path_length_cm,
                new_path_length_cm,
                optically_thin,
                waveunit,
                rescaled_units,
                extra,
                true_path_length,
            )
            apply_slit = apply_slit or slit_needed

        if "transmittance" in recompute:
            rescaled, rescaled_units, slit_needed = rescale_transmittance(
                spec,
                rescaled,
                initial,
                old_mole_fraction,
                new_mole_fraction,
                old_path_length_cm,
                new_path_length_cm,
                waveunit,
                rescaled_units,
                extra,
            )
            apply_slit = apply_slit or slit_needed

    if "xsection" in recompute:
        rescaled, rescaled_units = rescale_xsection(
            spec,
            rescaled,
            initial,
            old_mole_fraction,
            new_mole_fraction,
            old_path_length_cm,
            optically_thin,
            waveunit,
            rescaled_units,
            extra,
            true_path_length,
        )

    # delete former convoluted value if apply_slit will be used (just to be sure
    # we arent keeping a non rescaled value if something goes wrong)
    if apply_slit:
        for k in list(spec._q.keys()):
            if k != "wavespace" and k in CONVOLUTED_QUANTITIES:
                del spec._q[k]
                del spec.units[k]

    # Save (only) the ones that we want
    for q in wanted:
        if q in rescaled:
            spec._q[q] = rescaled[q]
            # # also add non-convoluted array until we apply slit
        elif apply_slit and q in CONVOLUTED_QUANTITIES and q + "_noslit" in rescaled:
            spec._q[q + "_noslit"] = rescaled[q + "_noslit"]

    # Update units
    spec.units.update(rescaled_units)

    # Reapply slit if needed
    # TODO: replace with directly convolving with slit stored in conditions
    # TODO: first, add an option to give arrays to apply_slit (see Issue #30)
    if apply_slit:
        if not (
            "slit_function" in spec.conditions
            and "slit_unit" in spec.conditions
            and "norm_by" in spec.conditions
        ):
            #            if 'transmittance' in extra and 'radiance' in extra:
            #                pass
            #            else:
            raise KeyError(
                "Slit is needed to recompute some quantities ({0}) ".format(wanted)
                + "but not all conditions are given among slit_function "
                + "({0}), slit_unit ({1}) and norm_by ({2})".format(
                    "slit_function" in spec.conditions,
                    "slit_unit" in spec.conditions,
                    "norm_by" in spec.conditions,
                )
            )
        else:
            slit_function = spec.conditions["slit_function"]
            slit_unit = spec.conditions["slit_unit"]
            norm_by = spec.conditions["norm_by"]
            try:
                shape = spec.conditions["slit_shape"]
            except KeyError:
                warn(
                    "Recomputing with slit_shape not given. Assuming a triangular slit.",
                    UserWarning,
                )
                shape = "triangular"
            spec.apply_slit(
                slit_function=slit_function,
                unit=slit_unit,
                shape=shape,
                norm_by=norm_by,
                verbose=verbose,
            )

        # Remove unwanted non-convoluted arrays
        for q in rescaled:
            if q in CONVOLUTED_QUANTITIES and q + "_noslit" not in wanted:
                del spec._q[q + "_noslit"]

    # Final checks

    # ... "everyone is here": check we didnt miss anyone
    rescaled_list = list(rescaled)
    # add the new quantities added by apply_slit
    rescaled_list = rescaled_list + [
        k for k in spec.get_vars() if k in CONVOLUTED_QUANTITIES
    ]
    for q in wanted:
        if not q in rescaled_list:
            raise AssertionError(
                "{0} could not be rescaled as wanted. ".format(q)
                + "The following properties were rescaled: {0}".format(rescaled_list)
            )
    # ... everyone was added in the Spectrum properly
    final_list = spec.get_vars()
    for q in wanted:
        if not q in final_list:
            raise AssertionError(
                "{0} is not in the final Spectrum. ".format(q)
                + "Rescaled spectrum contains: {0}".format(final_list)
            )
    # ... "everyone was rescaled": check we didnt scale only part of the spectrum
    for q in initial:
        if not q in rescaled_list:
            raise AssertionError(
                "{0} was initially in the Spectrum but was not ".format(q)
                + "rescaled. This can lead to error. Rescaled spectrum "
                + "contains: {0}".format(rescaled_list)
            )
    # ... all values have units
    for k in spec.get_vars():
        assert k in spec.units


def rescale_path_length(
    spec, new_path_length, old_path_length=None, inplace=False, force=False
):
    """Rescale spectrum to new path length. Starts from absorption coefficient
    and emission coefficient, and solves the RTE again for the new path length
    Convoluted values (with slit) are dropped in the process.

    Parameters
    ----------
    spec: Spectrum
    new_path_length: float (cm) or `~astropy.units.quantity.Quantity`
        new path length
    old_path_length: float (cm) or `~astropy.units.quantity.Quantity`, or None
        if None, current path length (conditions['path_length']) is used


    Other Parameters
    ----------------
    inplace: boolean
        if ``True``, modifies the Spectrum object directly. Else, returns
        a copy. Default ``False``.
    force: boolean
        if False, won't allow rescaling to 0 (not to loose information).
        Default ``False``

    Returns
    -------
    s_rescaled: Spectrum
        a rescaled Spectrum.
        if ``inplace=True``, then ``s`` has been rescaled already and
        ``s_rescaled`` is ``s``

    Notes
    -----
    Implementation:

        To deal with all the input cases, we first make a list of what has to
        be recomputed, and what has to be recalculated

    .. minigallery:: radis.Spectrum.rescale_path_length

    """

    if not inplace:
        spec = spec.copy()

    # Check inputs
    # ----------
    if old_path_length is not None:
        try:
            if spec.conditions["path_length"] != old_path_length:
                warn(
                    "path_length ({0}) doesnt match value given in conditions ({1})".format(
                        old_path_length, spec.conditions["path_length"]
                    )
                )
        except KeyError:  # path_length not defined
            pass
    else:
        try:
            old_path_length = spec.conditions["path_length"]
        except KeyError:
            raise KeyError(
                "path_length has to be defined in conditions (or use"
                + " `from_path_length`)"
            )

    if new_path_length < 0 and not force:
        raise ValueError("path_length cannot be negative")
    if new_path_length == 0 and not force:
        raise ValueError(
            "Rescaling to 0 will loose information. Choose force " "= True"
        )
    # Convert units
    new_path_length = convert_and_strip_units(new_path_length, u.cm)
    old_path_length = convert_and_strip_units(old_path_length, u.cm)

    for q in ["transmittance", "radiance"]:
        qns = q + "_noslit"
        qties = spec.get_vars()
        if q in qties and qns not in qties and not force:
            raise KeyError(
                "Cant rescale {0} if {1} not stored. ".format(q, qns)
                + " Use force=True to rescale anyway. {0}".format(q)
                + " will be deleted"
            )

    # Rescale
    assume_equilibrium = spec.conditions.get("thermal_equilibrium", False)
    _recalculate(
        spec,
        quantity="same",
        new_path_length_cm=new_path_length,
        old_path_length_cm=old_path_length,
        new_mole_fraction=1,
        old_mole_fraction=1,
        assume_equilibrium=assume_equilibrium,
    )

    # Update conditions
    spec.conditions["path_length"] = new_path_length

    return spec


def rescale_mole_fraction(
    spec,
    new_mole_fraction,
    old_mole_fraction=None,
    inplace=False,
    ignore_warnings=False,
    force=False,
    verbose=True,
):
    """Update spectrum with new molar fraction Convoluted values (with slit)
    are dropped in the process.

    Parameters
    ----------
    spec: Spectrum
    new_mole_fraction: float
        new mole fraction
    old_mole_fraction: float, or None
        if None, current mole fraction (conditions['mole_fraction']) is used


    Other Parameters
    ----------------
    inplace: boolean
        if ``True``, modifies the Spectrum object directly. Else, returns
        a copy. Default ``False``.
    force: boolean
        if False, won't allow rescaling to 0 (not to loose information).
        Default ``False``

    Returns
    -------
    s_rescaled: Spectrum
        a rescaled Spectrum.
        if ``inplace=True``, then ``s`` has been rescaled already and
        ``s_rescaled`` is ``s``

    Notes
    -----
    Implementation:

        similar to rescale_path_length() but we have to scale abscoeff & emisscoeff
        Note that this is valid only for small changes in mole fractions. Then,
        the change in line broadening becomes significant

    .. minigallery:: radis.Spectrum.rescale_mole_fraction

    """

    if not inplace:
        spec = spec.copy()

    # Check inputs
    # ---------
    # ... get old mole fraction, use existing one if not given
    if old_mole_fraction is not None:
        try:
            if (
                spec.conditions["mole_fraction"] != old_mole_fraction
                and not ignore_warnings
            ):
                warn(
                    "mole_fraction ({0}) doesnt match value given in conditions ({1})".format(
                        old_mole_fraction, spec.conditions["mole_fraction"]
                    )
                )
        except KeyError:  # mole fraction not defined
            pass

    else:
        try:
            old_mole_fraction = spec.conditions["mole_fraction"]
        except KeyError:
            raise KeyError(
                "mole_fraction has to be defined in conditions (or use"
                + " `from_mole_fraction`)"
            )

    # Add warning is large mole fraction rescale
    if np.abs(new_mole_fraction - old_mole_fraction) > 0.3:
        warn(
            "Large rescaling from {0} to {1}. ".format(
                old_mole_fraction, new_mole_fraction
            )
            + "There may be significant changes in pressure broadening coefficients."
            + "You should calculate a new spectrum instead."
        )

    if new_mole_fraction < 0 and not force:
        raise ValueError("mole_fraction cannot be negative")
    if new_mole_fraction == 0 and not force:
        raise ValueError(
            "Rescaling to 0 will loose information. Choose force " "= True"
        )
    if new_mole_fraction > 1 and not force:
        warn("rescaling to mole fraction > 1: {0}".format(new_mole_fraction))

    for q in ["transmittance", "radiance"]:
        qns = q + "_noslit"
        qties = spec.get_vars()
        if q in qties and qns not in qties and not force:
            raise KeyError(
                "Cant rescale {0} if {1} not stored.".format(q, qns)
                + "(you need to rescale before applying the slit again) "
                + " Use force=True to rescale anyway, but {0}".format(q)
                + " will be deleted"
            )

    # Get path length
    if "path_length" in list(spec.conditions.keys()):
        path_length = spec.conditions["path_length"]
        true_path_length = True
    else:
        path_length = (
            1  # path length not given, but we can still recompute a few spectral arrays
        )
        true_path_length = False

    # Rescale
    assume_equilibrium = spec.conditions.get("thermal_equilibrium", False)
    _recalculate(
        spec,
        quantity="same",
        new_path_length_cm=path_length,
        old_path_length_cm=path_length,
        new_mole_fraction=new_mole_fraction,
        old_mole_fraction=old_mole_fraction,
        true_path_length=true_path_length,
        assume_equilibrium=assume_equilibrium,
        verbose=verbose,
    )

    # Update conditions
    spec.conditions["mole_fraction"] = new_mole_fraction

    return spec


if __name__ == "__main__":

    from radis.test.spectrum.test_rescale import _run_all_tests

    print(("Test rescale.py: ", _run_all_tests(verbose=True)))

    def test_recompute_transmittance():

        import radis

        s = radis.test_spectrum()
        s.update("all")
        assert set(s.get_vars()) == set(
            [
                "abscoeff",
                "radiance_noslit",
                "absorbance",
                "transmittance_noslit",
                "emissivity_noslit",
                "emisscoeff",
                "xsection",
            ]
        )

        for k in list(s._q.keys()):
            if k == "wavespace":
                pass
            elif k == "abscoeff":
                pass
            else:
                del s._q[k]

        s.get("transmittance_noslit")
        print(s._q.keys())
        assert set(s.get_vars()) == set(("abscoeff", "transmittance_noslit"))
