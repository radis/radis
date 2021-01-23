# -*- coding: utf-8 -*-
"""Created on Tue Jul 18 17:47:48 2017.

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


# %% Dunham development

# ... vibrational term (Herzberg notation)


def Gv(v, we=0, wexe=0, weye=0, weze=0, weae=0, webe=0, gv=1):
    """Vibrational energy term Dunham development (order 5 in v) in Herzberg
    notation.

    .. math::

        G_v = w_e\\left(v+\\frac{g_v}{2}\\right) - w_ex_e\\left(v+\\frac{g_v}{2}\\right)^2 +
        w_ey_e\\left(v+\\frac{g_v}{2}\\right)^3  + w_ez_e\\left(v+\\frac{g_v}{2}\\right)^4 +
        w_ea_e\\left(v+\\frac{g_v}{2}\\right)^5

    Parameters
    ----------

    v: int
        vibrational quantum number
    we: float
            vibrational constant – first term (cm-1)
    ωexe: float
            vibrational constant – second term (cm-1)
    ωeye: float
            vibrational constant – third term (cm-1)
    ωeze: float
            vibrational constant – fourth term (cm-1)
    weae: float
       vibrational constant – fifth term (cm-1)
    webe:
        vibrational constant – sixth term (cm-1)
    gv: int
       degeneracy   (usually 1, but 2 for CO2-v2)

    Returns
    -------

    Gv: float
        Energy (cm-1)


    Notes
    -----

    Validity:

    For large vibrational levels Dunham's expansion is not valid. In other LBL codes,
    such as Specair, Morse Potential is used above a certain vibrational level.
    See :func:`~radis.phys.morse.morse_increment` for an implementation in RADIS.

    .. warning::

        Although by convention prefactors share the same name throughout most of the litterature,
        signs can be different depending on the article. Make sure your Dunham expansion
        has the same signs as the one we use here! In particular, this function
        uses a ``-wexe`` sign as in the usual Herzberg notation. Ex:

        - Mantz and Maillard 1975 (CO X)  uses opposite signs for **weze** and **webe**

    References
    ----------

    "Optical Diagnostics and Radiative Emission of Air Plasmas", C. Laux, 1993, p93

    See Also
    --------

    :py:func:`~radis.levels.dunham.Fv`

    List of :py:data:`~radis.db.conventions.herzberg_coefficients`

    Get Herzberg coefficients for a molecule: :py:func:`~radis.db.utils.get_herzberg_coefficients`

    Use Dunham coefficients instead: :py:func:`~radis.levels.dunham.EvJ`
    """

    return (
        we * (v + gv / 2)
        - wexe * (v + gv / 2) ** 2
        + weye * (v + gv / 2) ** 3
        + weze * (v + gv / 2) ** 4
        + weae * (v + gv / 2) ** 5
        + webe * (v + gv / 2) ** 6
    )  # + ..
    #                :                                :
    #                :                                :
    # In Mantz and Maillard 1975                      - weze*(v+gv/2)**4
    # Mantz 1975     - webe*(v+gv/2)**6


# ... rotational term (Herzberg notation)


def Fv(
    v,
    J,
    Be=0,
    De=0,
    alpha_e=0,
    beta_e=0,
    gamma_e=0,
    delta_e=0,
    epsilon_e=0,
    pi_e=0,
    He=0,
    eta_e=0,
    gv=1,
):
    """Rotational energy term Dunham development (order 4 in J) in Herzberg
    notation.

    .. math::
        B_{v}=B_{e}-\\alpha_{e}\\left(v+\\frac{g_{v}}{2}\\right)+\\gamma_{e}
        \\left(v+\\frac{g_{v}}{2}\\right)^{2}+\\delta_{e}\\left(v+\\frac{g_{v}}{2}
        \\right)^{3}+\\epsilon_{e}\\left(v+\\frac{g_{v}}{2}\\right)^{4}

        D_{v}=D_{e}+\\beta_{e}\\left(v+\\frac{g_{v}}{2}\\right)+\\pi_{e}
        \\left(v+\\frac{g_{v}}{2}\\right)^{2}

        H_{v}=H_{e}-\\eta_{e}\\left(v+\\frac{g_{v}}{2}\\right)

    *generated from the Python formula with* :py:func:`~pytexit.pytexit.py2tex`


    Parameters
    ----------

    v: int
        vibrational quantum number
    J: int
        rotational quantum number
    Be: float
            rotational constant in equilibrium position (cm-1)
    De: float
            centrifugal distortion constant (cm-1)
    alpha_e: float
            rotational constant – first term (cm-1)
    beta_e: float
            rotational constant – first term, centrifugal force (cm-1)
    gamma_e: float
            rotation-vibration interaction constant (cm-1)
    delta_e: float
            (cm-1)
    epsilon_e: float
       (cm-1)
    pi_e: float
        (cm-1)
    He: float
        third order correction factor (cm-1)
    eta_e: float
        (float)
    gv: int
       degeneracy   (usually 1, but 2 for CO2-v2)

    Returns
    -------

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

    See Also
    --------

    :py:func:`~radis.levels.dunham.Gv`

    List of :py:data:`~radis.db.conventions.herzberg_coefficients`

    Get Herzberg coefficients for a molecule: :py:func:`~radis.db.utils.get_herzberg_coefficients`

    Use Dunham coefficients instead: :py:func:`~radis.levels.dunham.EvJ`
    """

    B_e, D_e, H_e, g_v = Be, De, He, gv

    # ... Note: formula added in docstring with pytexit.py2tex()
    # + ...
    B_v = (
        B_e
        - alpha_e * (v + g_v / 2)
        + gamma_e * (v + g_v / 2) ** 2
        + delta_e * (v + g_v / 2) ** 3
        + epsilon_e * (v + g_v / 2) ** 4
    )
    D_v = D_e + beta_e * (v + g_v / 2) + pi_e * (v + g_v / 2) ** 2  # + ...
    H_v = H_e - eta_e * (v + g_v / 2)
    #       :                                        :
    #       :                                        :
    # In Mantz and Maillard 1975                     - delta_e*(v+gv/2)**3
    #       - beta_e*(v+gv/2)      in Mantz and Maillard 1975

    return B_v * J * (J + 1) - D_v * (J * (J + 1)) ** 2 + H_v * (J * (J + 1)) ** 3


# ... general term

# def EvJ(v, J, **Ykl_dict):
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
#        from radis.db.utils import get_dunham_coefficients
#        from radis.levels.dunham import EvJ
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
    """Calculates rovibrational energy reading from Dunham coefficients in Ykl
    notation, for diatomic molecules.

    Parameters
    ----------

    Ykl: dict
        an arbitrary dictionary of Ykl coefficients
        accepted formats: Y01

    Ykl are parsed and assigned the correct energy

    Examples
    --------

    Read directly from a .json file::

        from radis.db.utils import get_dunham_coefficients
        from radis.levels.dunham import EvJ
        dunham_coeffs = get_dunham_coefficients('CO', 1, 'X1SIG+')

        # Now calculate energy
        EvJ(v=0, J=0, **dunham_coeffs)

    See Also
    --------

    Get Dunham coefficients for a molecule: :py:func:`~radis.db.utils.get_dunham_coefficients`

    Use Herzberg coefficients instead: :py:func:`~radis.levels.dunham.Gv`,
    :py:func:`~radis.levels.dunham.Fv`
    """

    E = 0

    for ykl_name, Ykl in Ykl_dict.items():
        k = int(ykl_name[1])
        l = int(ykl_name[2])
        E += Ykl * (v + 0.5) ** k * (J * (J + 1)) ** l

    return E


if __name__ == "__main__":

    from radis.test.phys.test_dunham import _run_all_tests

    print("Testing Dunham.py: ", _run_all_tests())
