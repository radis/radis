# -*- coding: utf-8 -*-
"""Created on Mon Nov 13 21:58:47 2017.

@author: erwan

Formula to calculate most spectral quantities

Are stored here all equations that need to be used by both the factory (to
generate the spectrum) or the spectrum class (to recompute some quantities
in post processing)
"""

from radis.phys.blackbody import planck
from radis.phys.convert import cm2nm

# ------- Get radiance from emissivity


def calc_radiance(wavenumber, emissivity, Tgas, unit="mW/sr/cm2/nm"):
    """Derive radiance (:math:`mW/cm^2/sr/nm`) from the emissivity.

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

    radiance = emissivity * planck(cm2nm(wavenumber), Tgas, unit=unit)

    return radiance
