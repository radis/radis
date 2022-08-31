# -*- coding: utf-8 -*-
"""Created on Mon Nov 13 21:58:47 2017.

@author: erwan

Formula to calculate most spectral quantities

Are stored here all equations that need to be used by both the factory (to
generate the spectrum) or the spectrum class (to recompute some quantities
in post processing)
"""

from radis.phys.blackbody import planck, planck_wn
from radis.phys.constants import k_b
from radis.phys.convert import cm2nm
from radis.phys.units import is_homogeneous

# ------- Get radiance from emissivity


def calc_radiance(wavenumber, emissivity, Tgas, unit="mW/sr/cm2/nm"):
    r"""Derive radiance (:math:`mW/cm^2/sr/nm`) from the emissivity.

    .. math::
        I(\lambda) = \epsilon(\lambda) \cdot L_0(\lambda, T)

    Returns
    -------

    radiance: array in :math:`mW/sr/cm2/nm`
        unless ``unit`` is different

    See Also
    --------

    :py:func:`~radis.phys.blackbody.planck`
    """

    if is_homogeneous(unit, "mW/sr/cm2/nm"):
        radiance = emissivity * planck(cm2nm(wavenumber), Tgas, unit=unit)
    elif is_homogeneous(unit, "mW/sr/cm2/cm-1"):
        radiance = emissivity * planck_wn(wavenumber, Tgas, unit=unit)

    return radiance


def abscoeff2xsection(abscoeff_cm1, Tgas_K, pressure_Pa, mole_fraction):
    r"""Convert Absorption coefficient (cm-1) to Cross-section  (cm2)

    Parameters
    ----------
    abscoeff_cm : cm-1
    Tgas_K : K
    pressure_Pa : Pa
    mole_fraction : 0-1

    Returns
    -------
    xsection: cm2
    """

    return abscoeff_cm1 * (k_b * Tgas_K / mole_fraction / pressure_Pa) * 1e6  # cm2


if __name__ == "__main__":
    from radis.test.spectrum.test_equations import test_equations

    test_equations()
