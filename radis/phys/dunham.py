# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 17:47:48 2017

@author: erwan

Dunham development for diatomic molecules energies

Warning
-------

Although by convention prefactors share the same name throughout most of the litterature,
signs can be different depending on the article. Make sure your Dunham expansion
has the same signs as the one we use here!

Reference
------------
"Optical Diagnostics and Radiative Emission of Air Plasmas", C. Laux, 1993, p93
Mantz et al 1975, "Ground state molecular constants of 12C16O" 
    
"""

from __future__ import division, absolute_import, print_function
from __future__ import unicode_literals  # α, β, γ, δ, ϵ
import re

# %% Dunham development

# ... vibrational term (Herzberg notation)

def Gv(v, we, wexe, weye, weze, weae=0, webe=0, gv=1):
    ''' Vibrational energy term 
    Dunham development (order 5 in v) in Herzberg notation
    
    >>> Gv = we*(v+gv/2) - wexe*(v+gv/2)^2 + weye*(v+gv/2)^3 + weze*(v+gv/2)^4 + 
    >>> weae*(v+gv/2)^5
        
    Parameters
    ----------
    
    v        vibrational quantum number
    we 	    vibrational constant – first term (cm-1)
    ωexe 	    vibrational constant – second term (cm-1)
    ωeye 	    vibrational constant – third term (cm-1)
    ωeze 	    vibrational constant – fourth term (cm-1)
    weae     vibrational constant – fifth term (cm-1)
    webe     vibrational constant – sixth term (cm-1)
    gv       degeneracy   (usually 1, but 2 for CO2-v2) 
    
    Returns
    -------
    
    Gv: float
        Energy (cm-1)


    Notes
    -----
    
    Validity:

    For large vibrational levels Dunham's expansion is not valid. In Specair
    Morse Potential is used above a certain vibrational level

    .. warning::

        Although by convention prefactors share the same name throughout most of the litterature,
        signs can be different depending on the article. Make sure your Dunham expansion
        has the same signs as the one we use here! Ex:
    
        - Mantz and Maillard 1975 (CO)  uses opposite signs for **weze** and **webe**

    References
    ----------

    "Optical Diagnostics and Radiative Emission of Air Plasmas", C. Laux, 1993, p93

    '''

    return we*(v+gv/2) - wexe*(v+gv/2)**2 + weye*(v+gv/2)**3 + weze*(v+gv/2)**4 + \
           weae*(v+gv/2)**5 + webe*(v+gv/2)**6 # + .. 
           #                :                                :
           #                :                                :
           # In Mantz and Maillard 1975                      - weze*(v+gv/2)**4
           # Mantz 1975     - webe*(v+gv/2)**6
           
# ... rotational term (Herzberg notation)

def Fv(v, J, Be, De, alpha_e, beta_e, gamma_e=0, delta_e=0, epsilon_e=0,
       pi_e=0, He=0, eta_e=0, gv=1):
    '''Rotational energy term
    Dunham development (order 4 in J) in Herzberg notation
    
    .. math::
        B_{v}=B_{e}-\\alpha_{e}\\left(v+\\frac{g_{v}}{2}\\right)+\\gamma_{e}
        \\left(v+\\frac{g_{v}}{2}\\right)^{2}+\\delta_{e}\\left(v+\\frac{g_{v}}{2}
        \\right)^{3}+\\epsilon_{e}\\left(v+\\frac{g_{v}}{2}\\right)^{4}

        D_{v}=D_{e}+\\beta_{e}\\left(v+\\frac{g_{v}}{2}\\right)+\\pi_{e}
        \\left(v+\\frac{g_{v}}{2}\\right)^{2}

        H_{v}=H_{e}-\\eta_{e}\\left(v+\\frac{g_{v}}{2}\\right)

    *generated from code with pytexit*


    Parameters
    ----------

    v        vibrational quantum number
    J        rotational quantum number
    Be 	    rotational constant in equilibrium position (cm-1)
    De 	    centrifugal distortion constant (cm-1)
    alpha_e 	    rotational constant – first term (cm-1)
    beta_e 	    rotational constant – first term, centrifugal force (cm-1)
    gamma_e 	    rotation-vibration interaction constant (cm-1)
    delta_e 	    (cm-1)
    epsilon_e       (cm-1)
    pi_e
    He      third order correction factor
    eta_e      

    gv       degeneracy   (usually 1, but 2 for CO2-v2) 

    Returns
    ------
    
    Fv: float
        Energy (cm-1)
    
    Notes
    -----
    
    Validity:

    For large vibrational levels Dunham's expansion is not valid. In RADIS a 
    Morse Potential can be used above a certain vibrational level

    .. warning::
    
        Although by convention prefactors share the same name throughout most of the litterature,
        signs can be different depending on the article. Make sure your Dunham expansion
        has the same signs as the one we use here! Ex:
    
        - Mantz and Maillard 1975  (CO)  uses opposite signs for **delta_e** and **beta_e**

    References
    ----------

    "Optical Diagnostics and Radiative Emission of Air Plasmas", C. Laux, 1993, p93

    '''

    B_e, D_e, H_e, g_v = Be, De, He, gv

    # ... Note: formula added in docstring with pytexit.py2tex()
    # + ...
    B_v = B_e - alpha_e*(v+g_v/2) + gamma_e*(v+g_v/2)**2 + \
        delta_e*(v+g_v/2)**3 + epsilon_e*(v+g_v/2)**4
    D_v = D_e + beta_e*(v+g_v/2) + pi_e*(v+g_v/2)**2  # + ...
    H_v = H_e - eta_e*(v+g_v/2)
    #       :                                        :
    #       :                                        :
    # In Mantz and Maillard 1975                     - delta_e*(v+gv/2)**3
    #       - beta_e*(v+gv/2)      in Mantz and Maillard 1975

    
    return B_v*J*(J+1) - D_v*(J*(J+1))**2 + H_v*(J*(J+1))**3


# ... general term
    
#def EvJ(v, J, **Ykl_dict):
#    ''' Calculates rovibrational energy reading from Dunham coefficients in 
#    Ykl notation
#    
#    Parameters
#    ----------
#    
#    Ykl: dict
#        an arbitrary dictionary of Ykl coefficients
#        accepted formats: Y01, Y01_cm-1
#        
#    Ykl are parsed and assigned the correct energy
#    
#    Notes
#    -----
#    
#    Because reading and parsing the Ykl dict is required, this is expected to 
#    be slower that a hardcoded implementation. However, it is also much 
#    more flexible. 
#
#    
#    Examples
#    --------
#    
#    Read directly from a .json file::
#            
#        from neq.db.utils import get_dunham_coefficients
#        from neq.phys.dunham import EvJ
#        dunham_coeffs = get_dunham_coefficients('CO', 1, 'X1SIG+')
#        
#        # Now calculate energy
#        EvJ(v=0, J=0, **dunham_coeffs)
#    
#    '''
#    
#    E = 0
#    
#    Ykl_format = '^Y(?P<k>[0-9])(?P<l>[0-9])(_cm-1)?$'
#    # ... Ykl with k, l ints  and followed by nothing or _cm-1 
#    # ... accepted formats: Y01, Y01_cm-1
#    reg = re.compile(Ykl_format)
#    
#    for ykl_name, Ykl in Ykl_dict.items():
#        match = reg.search(ykl_name)
#        if match is None:
#            raise ValueError('Key {0} does not have the expected format {1}'.format(
#                    ykl_name, Ykl_format))
#        res = match.groupdict()
#        k = int(res['k'])
#        l = int(res['l'])
#        E += Ykl * (v+0.5)**k * (J*(J+1))**l
#    
#    return E

# New version (less versatile, but less Python and faster)
# only allowed formats: Y01
def EvJ(v, J, **Ykl_dict):
    ''' Calculates rovibrational energy reading from Dunham coefficients in 
    Ykl notation
    
    Parameters
    ----------
    
    Ykl: dict
        an arbitrary dictionary of Ykl coefficients
        accepted formats: Y01
        
    Ykl are parsed and assigned the correct energy
    
    Examples
    --------
    
    Read directly from a .json file::
            
        from neq.db.utils import get_dunham_coefficients
        from neq.phys.dunham import EvJ
        dunham_coeffs = get_dunham_coefficients('CO', 1, 'X1SIG+')
        
        # Now calculate energy
        EvJ(v=0, J=0, **dunham_coeffs)
    
    '''
    
    E = 0
    
    for ykl_name, Ykl in Ykl_dict.items():
        k = int(ykl_name[1])
        l = int(ykl_name[2])
        E += Ykl * (v+0.5)**k * (J*(J+1))**l
    
    return E

if __name__ == '__main__':
    
    from neq.test.phys.test_dunham import _run_all_tests
    print('Testing Dunham.py: ', _run_all_tests())
    