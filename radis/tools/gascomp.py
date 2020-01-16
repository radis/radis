# -*- coding: utf-8 -*-
"""
Determine Gas composition using [CANTERA]_

References
----------

.. [CANTERA]  D. G. Goodwin, H. K. Moffat, R. L. Speth, Cantera: An Object-oriented Software
              Toolkit for Chemical Kinetics, Thermodynamics, and Transport Processes,
              http://www.cantera.org, doi:10.5281/zenodo.170284, 2017.

-------------------------------------------------------------------------------

"""

from __future__ import absolute_import
from radis.misc.utils import NotInstalled

try:
    import cantera as ct
except:
    ct = NotInstalled(
        "cantera",
        "Cantera is needed to calculate equilibrium mole fractions"
        + ". Install with  `pip install cantera`",
    )


def get_eq_mole_fraction(initial_mixture, T_K, p_Pa):
    """ Returns mole fraction at temperature T, using the 
    [CANTERA]_ :py:meth:`~cantera.ThermoPhase.equilibrate` function. 
    
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

    # %% Init Cantera
    g = ct.Solution("gri30.xml")
    g.TPX = T_K, p_Pa, initial_mixture

    # Minimize Gibbs:
    g.equilibrate("TP")

    # Returns
    return g.mole_fraction_dict()
