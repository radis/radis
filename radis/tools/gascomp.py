# -*- coding: utf-8 -*-
"""
Determine gas mixture composition under chemical equilibrium using CANTERA.

-------------------------------------------------------------------------------

"""


def get_eq_mole_fraction(initial_mixture, T_K, p_Pa):
    """Calculates chemical equilibrium mole fraction at temperature T, using
    the CANTERA :py:meth:`~cantera.ThermoPhase.equilibrate` function.

    The calculation uses the default GRI3.0 mechanism, which was
    designed to model natural gas combustion, including NO formation
    and reburn chemistry. See `GRI 3.0 <combustion.berkeley.edu/gri-mech/version30/text30.html>`__.

    When using, cite the [CANTERA]_ package.

    Parameters
    ----------

    initial_mixture: str
        Gas composition. Example::

             'N2:0.79, O2:0.21, CO2:363e-6'

        Or::

             'CO2:1'

    T_K: float (K)
        temperature (Kelvin) to calculate equilibrium

    P_Pa: float (Pa)
        temperature (Pascal) to calculate equilibrium

    Examples
    --------

    Calculate equilibrium mixture of CO2 at 2000 K, 1 atm::

        get_eq_mole_fraction('CO2:1', 2000, 101325)

        >>> {'C': 1.7833953335281855e-19,
            'CO': 0.01495998583472384,
            'CO2': 0.9775311634424326,
            'O': 5.7715610124613225e-05,
            'O2': 0.007451135112719029}

    References
    ----------

    [CANTERA]_
    """

    try:
        import cantera as ct
    except ImportError as err:
        raise ImportError(
            "Cantera is needed to calculate equilibrium mole fractions"
            + ". Install with  `pip install cantera` or (better) `conda install -c cantera cantera`",
        ) from err

    # %% Init Cantera
    g = ct.Solution("gri30.xml")
    g.TPX = T_K, p_Pa, initial_mixture

    # Minimize Gibbs:
    g.equilibrate("TP")

    # Returns
    return g.mole_fraction_dict()
