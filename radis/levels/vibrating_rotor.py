# -*- coding: utf-8 -*-
"""Created on Tue Jul 18 17:59:45 2017.

@author: erwan


Rovibrational Energy for Vibrating Rotor

Ex: used for CO2 X(1Σ+) state

References
----------





-----------------------------------------------------------------------------
"""


from radis.db.conventions import (
    herzberg_coefficients_rot,
    herzberg_coefficients_rovib,
    herzberg_coefficients_vib,
)
from radis.levels.dunham import Fv, Gv


def EvJ_uncoupled_vibrating_rotor(
    v1, v2, l2, v3, J, coeff_dict, gv1=1, gv2=1, gv3=1, remove_ZPE=True
):
    """Rovibrational energy of an uncoupled vibrating rotor.

    Parameters
    ----------

    v1: int
        vibrational state

    v2: int
        vibrational state

    l2: int
        vibrational state

    v3: int
        vibrational state

    J: int
        rotational state

    coeff_dict: dict
        dictionary of :py:data:`~radis.db.conventions.herzberg_coefficients`,
        which can be numbered for the different vibration modes. Example::

            {'we1': 1333.93,
             'we2': 667.47,
             'we3': 2349.16,
             'wexe1': 2.93,
             'wexe2': -0.38,
             'wexe3': 12.47,
             'Be': 0.39022,
             'De': 1.333e-07,
             'He': 9e-15}

    gv1, gv2, gv3: int
        degeneracies of each vibrational mode. Default ``1, 1, 1``::

            1,2,1 for CO2

    remove_ZPE: boolean
        if ``True``, removes energy of ground state vibrational level (zero-point-energy)

    Returns
    -------

    E: float
        energy of state in cm-1

    References
    ----------

    Klarenaar et al, "Time evolution of vibrational temperatures in a CO 2 glow
    discharge measured with infrared absorption spectroscopy", doi 10.1088/1361-6595/aa902e,
    and the references there in.

    See Also
    --------

    :py:func:`~radis.levels.energies_co2.EvJ_co2`
    """

    coeff_dict = coeff_dict.copy()

    coeffs_vib = {"1": {}, "2": {}, "3": {}}
    coeffs_rot = {}

    # Split coeffs between the different vibration and rotation modes:
    for k, v in coeff_dict.items():
        if k in herzberg_coefficients_rot:
            coeffs_rot[k] = v
        elif k[:-1] in herzberg_coefficients_vib:
            vib_mode = k[-1]
            coeffs_vib[vib_mode][k[:-1]] = v
        elif k in herzberg_coefficients_rovib:
            raise NotImplementedError(
                "Mixed term {0} not implemented for uncoupled rovib model".format(k)
            )
        else:
            raise KeyError("Unexpected coefficient: {0}".format(k))

    coeffs_vib1 = coeffs_vib["1"]
    coeffs_vib2 = coeffs_vib["2"]
    coeffs_vib3 = coeffs_vib["3"]

    # Get rotational coeffs:
    coeffs_rot = {
        k: v for (k, v) in coeff_dict.items() if k in herzberg_coefficients_rot
    }

    # Energies
    G1 = Gv(v1, gv=gv1, **coeffs_vib1)
    G2 = Gv(v2, gv=gv2, **coeffs_vib2)
    G3 = Gv(v3, gv=gv3, **coeffs_vib3)
    G = G1 + G2 + G3
    F = Fv(0, J, **coeffs_rot)  # Uncoupled model: v dependance ignored

    if remove_ZPE:
        ZPE = Gv(0, **coeffs_vib1) + Gv(0, gv=2, **coeffs_vib2) + Gv(0, **coeffs_vib3)
    else:
        ZPE = 0

    return G + F - ZPE  # cm-1


def EvJah_uncoupled_vibrating_rotor(v1, v2, l2, v3, J, coeff_dict, remove_ZPE=True):
    """Return rovibrationalenergies for nu_1, nu_2, nu_3 vibrational modes Each
    energy is a tuple (E_harmonic, E_nonharmonic) to be used for instance in a
    Treanor distribution.

    Defined for CO2 only. See:  Klarenaar et al. [1]_

    Energy of CO2 X(1Σg+) state with an independant harmonic oscillator
    approximation  (CO2-626 = first isotope)

    Note that this is less accurate than using CDSD energies in the (p, J, c, n)
    general denomination

    Parameters
    ----------

    v1: int
        vibrational state

    v2: int
        vibrational state

    l2: int
        vibrational state

    v3: int
        vibrational state

    J: int
        rotational state

    remove_ZPE: boolean
        if ``True``, removes energy of ground state vibrational level (zero-point-energy)

    Returns
    -------

    E: float
        energy of state in cm-1

    References
    ----------

    .. [1] Eq.(10) in Klarenaar et al, "Time evolution of vibrational temperatures in a CO2 glow
        discharge measured with infrared absorption spectroscopy", doi 10.1088/1361-6595/aa902e,
        and the references there in.
    """

    coeff_dict = coeff_dict.copy()

    we1 = coeff_dict.pop("we1")
    we2 = coeff_dict.pop("we2")
    we3 = coeff_dict.pop("we3")

    wexe1 = coeff_dict.pop("wexe1")
    wexe2 = coeff_dict.pop("wexe2")
    wexe3 = coeff_dict.pop("wexe3")

    Be = coeff_dict.pop("Be")
    De = coeff_dict.pop("De")
    He = coeff_dict.pop("He")

    if len(coeff_dict) > 0:
        raise NotImplementedError(
            "Harmonic/anharmonic energy split not defined "
            + "with the following spectroscopic constants: {0}".format(
                list(coeff_dict.keys())
            )
        )

    # Vibrational energy
    G1_h = we1 * v1
    G1_a = -wexe1 * v1 * (v1 - 1)
    # G1 = G1_h + G1_a

    G2_h = we2 * v2
    G2_a = -wexe2 * v2 * (v2 - 1)
    # G2 = G2_h + G2_a

    G3_h = we3 * v3
    G3_a = -wexe3 * v3 * (v3 - 1)
    # G3 = G3_h + G3_a

    # G = G1 + G2 + G3

    # Rotational energy
    Bv = Be  # no coupling terms
    Dv = De
    Hv = He
    F = Bv * J * (J + 1) - Dv * J ** 2 * (J + 1) ** 2 + Hv * J ** 3 * (J + 1) ** 3

    if remove_ZPE:
        # ZPE = 0  # by construction here
        pass
    else:
        raise NotImplementedError

    return (G1_h, G1_a), (G2_h, G2_a), (G3_h, G3_a), F


#    return Te + G + F - ZPE   # cm-1


# %% Test


def _test(*args, **kwargs):
    #    E_CO2X = E_CO2_626_X1SIGg
    from radis.db.utils import get_herzberg_coefficients

    coeffs = get_herzberg_coefficients("CO2", 1, "X1SIGu+")
    E_CO2X = lambda v1, v2, l2, v3, J: EvJ_uncoupled_vibrating_rotor(
        v1, v2, l2, v3, J, coeff_dict=coeffs, gv1=1, gv2=2, gv3=1
    )
    print(
        "CO2((0,0,0,0),J=0) -> CO2((0,0,0,1),J=0) energy: ",
        E_CO2X(0, 0, 0, 1, 0) - E_CO2X(0, 0, 0, 0, 0),
        "cm-1",
    )


if __name__ == "__main__":

    _test()
