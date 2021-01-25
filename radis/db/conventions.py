# -*- coding: utf-8 -*-
"""Tools to convert spectroscopic constants with the conventions of Hubert and
Herzberg (ω, Be, α, β, etc.) to Dunham parameters (Ykl)

References
----------

https://en.wikipedia.org/wiki/Dunham_expansion
"""

dunham2herzberg = {
    "Y01": (1, "Be"),
    "Y02": (-1, "De"),
    "Y03": (1, "He"),
    "Y04": (1, "Le"),
    "Y10": (1, "we"),
    "Y11": (-1, "alpha_e"),
    "Y12": (-1, "beta_e"),
    "Y20": (-1, "wexe"),
    "Y21": (-1, "gamma_e"),
    "Y30": (1, "weye"),
    "Y40": (1, "weze"),
}
""" dict: {Yij: (sign, coeff) }
conversion of Dunham spectroscopic coefficients to Herzberg convention """

## Invert the dictionary.
herzberg2dunham = {v: (sign, k) for k, (sign, v) in dunham2herzberg.items()}
""" dict: {Yij: (sign, coeff)
conversion of Herberg convention to Dunham spectroscopic coefficients

See Also
--------

:py:func:`~radis.db.utils.get_dunham_coefficients`,
:py:func:`~radis.db.utils.get_herzberg_coefficients`
"""

# Name of all Herzberg coefficients
# herzberg_coefficients_elec = ['Te']
#'''list: Herzberg coefficients used to calculate pure electronic energies'''
herzberg_coefficients_vib = ["we", "wexe", "weye", "weze", "weae", "webe"]
"""list: Herzberg coefficients used to calculate pure vibrational energies"""
herzberg_coefficients_rot = ["Be", "De", "He", "Le"]
"""list: Herzberg coefficients used to calculate pure rotational energies"""
herzberg_coefficients_rovib = [
    "alpha_e",
    "beta_e",
    "gamma_e",
    "delta_e",
    "eta_e",
    "pi_e",
]
"""list: Herzberg coefficients used in vibrational-rotational coupling terms for
rovibrational energies"""

herzberg_coefficients = (
    herzberg_coefficients_vib + herzberg_coefficients_rot + herzberg_coefficients_rovib
)
"""list: Herzberg coefficients used to calculate rovibrational energies in
:func:`~radis.levels.dunham.Fv`, :func:`~radis.levels.dunham.Gv`

See Also
--------

:py:func:`~radis.db.utils.get_herzberg_coefficients`

"""


# Sanity check (test coefficients defined in conversion dictionaries are valid
# Dunham and Herzberg coefficients)
for k in dunham2herzberg:
    assert k.startswith("Y")
for k in herzberg2dunham:
    assert k in herzberg_coefficients


def get_convention(coefficients):
    """Returns if we're using the Herzberg or Dunham convention for
    spectrosopic coefficients, and returns the associated coefficients.

    Parameters
    ----------

    coefficients:
        list

    Returns
    -------

    convention: str
        'dunham' or 'herzberg'

    coefficients: dict
        list of coefficient names
    """

    from radis.db.utils import ignore_trailing_number
    from radis.misc.basics import partition

    assert len(coefficients) > 0

    herzberg_coeffs, non_herzberg_coeffs = partition(
        lambda x: ignore_trailing_number(x) in herzberg_coefficients, coefficients
    )
    for k in non_herzberg_coeffs:
        if not k.startswith("Y"):
            raise ValueError("Unexpected Spectroscopic coefficient: {0}".format(k))
    dunham_coeffs = non_herzberg_coeffs

    if len(dunham_coeffs) > 0 and len(herzberg_coeffs) > 0:
        raise ValueError(
            "Both Dunham ({0}) and Herzberg ({1}) conventions used".format(
                dunham_coeffs, herzberg_coeffs
            )
            + ". Choose one only"
        )

    if len(dunham_coeffs) > 0:
        return "dunham"
    else:
        return "herzberg"


if __name__ == "__main__":

    assert get_convention(["wexe"]) == "herzberg"
    assert get_convention(["wexe1"]) == "herzberg"
    assert get_convention(["Y01", "Y11"]) == "dunham"
    import pytest

    with pytest.raises(ValueError):
        # mixed conventions: should fail
        get_convention(["wexe", "Y01", "Y11"])
