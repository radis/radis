# -*- coding: utf-8 -*-
"""Unique file where to define names.

Summary
-------

Standardize names for:

- vibrational levels
- ...



Routine Listing
---------------

vib_lvl_name_hitran_class1,
vib_lvl_name_hitran_class5,
vib_lvl_name_hitran_class5_short,
vib_lvl_name_cdsd_p,
vib_lvl_name_cdsd_pc,
vib_lvl_name_cdsd_pcN


----------
"""

# %% Vibrational levels names
# for all conventions (HITRAN, CDSD, etc.)

# Note for developers:
# ... dont use '..'.format() syntax here else it will break compatibility
# ... when working with Pandas columns.as_str() as input

# HITRAN ("spectroscopic") convention


def vib_lvl_name_hitran_class1(v1):
    """Write vibrational level for a HITRAN class 1 molecule: see
    :data:`~radis.io.hitran.HITRAN_CLASS1`
    """
    (v1,) = _format_str(v1)
    return "(" + v1 + ")"


def vib_lvl_name_hitran_class5(v1, v2, l2, v3):
    """Write vibrational level for a HITRAN class 5 molecule: see
    :data:`~radis.io.hitran.HITRAN_CLASS5`
    """
    v1, v2, l2, v3 = _format_str(v1, v2, l2, v3)
    return "(" + v1 + "," + v2 + "," + l2 + "," + v3 + ")"


def vib_lvl_name_hitran_class5_short(v1, v2, l2, v3):
    """Write vibrational level for a HITRAN class 5 molecule: see
    :data:`~radis.io.hitran.HITRAN_CLASS5`
    """
    v1, v2, l2, v3 = _format_str(v1, v2, l2, v3)
    return v1 + v2 + "`" + l2 + "`" + v3


# CDSD vibrational conventions (for CO2 only)


def vib_lvl_name_cdsd_p(p):
    """Write vibrational level with CDSD format
    Convention: we use (p,c)   (polyad, wang)
    """
    (p,) = _format_str(p)  # , is important
    return "(" + p + ")"


def vib_lvl_name_cdsd_pc(p, c):
    """Write vibrational level with CDSD format
    Convention: we use (p,c)   (polyad, wang)
    """
    p, c = _format_str(p, c)
    return "(" + p + "," + c + ")"


def vib_lvl_name_cdsd_pcN(p, c, N):
    """Write vibrational level with CDSD format
    Convention: we use (p,c,N)   (polyad, wang, rank number)

    Notes
    -----

    the N quantum number carries little physical sense: it's just a ranking number
    within a (p,c,J) group. But, it makes a (p,c,N), (J) level unique.
    """
    p, c, N = _format_str(p, c, N)
    return "(" + p + "," + c + "," + N + ")"


def vib_lvl_name_cdsd_pcJN(p, c, J, N):
    """Write vibrational level with CDSD format
    Convention: vibrational energy defined uniquely for all levels:
        (p,c,J,N)   (polyad, wang, rotational number, rank number)

    Notes
    -----

    the N quantum number carries little physical sense: it's just a ranking number
    within a (p,c,J) group. But, it makes a (p,c,N), (J) level unique.
    """
    p, c, J, N = _format_str(p, c, J, N)
    return "(" + p + "," + c + "," + J + "," + N + ")"


# %% Utils


from pandas import Series


def _format_str(*var):
    """Convert variables ``var`` to str.

    Uses ``.astype(str)`` if they are Pandas series
    """

    out = []
    for v in var:
        if isinstance(v, Series):
            out.append(v.astype(str))
        else:
            out.append(str(v))
    return out


if __name__ == "__main__":

    from radis.test.lbl.test_labels import test_vibrational_levels_labelling

    test_vibrational_levels_labelling()
