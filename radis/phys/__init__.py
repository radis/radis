# -*- coding: utf-8 -*-
"""Physical constants and conversion.
"""

from .blackbody import planck, planck_wn, sPlanck
from .convert import (
    J2K,
    K2J,
    J2cm,
    J2eV,
    K2eV,
    atm2bar,
    atm2torr,
    bar2atm,
    bar2torr,
    cm2eV,
    cm2hz,
    cm2J,
    cm2nm,
    dcm2dnm,
    dnm2dcm,
    eV2cm,
    eV2J,
    eV2K,
    eV2nm,
    hz2cm,
    nm2cm,
    nm2eV,
    torr2atm,
    torr2bar,
)
from .units import (
    conv2,
    convert_emi2cm,
    convert_emi2nm,
    convert_rad2cm,
    convert_rad2nm,
    is_homogeneous,
)
from .units_astropy import convert_and_strip_units

__all__ = ["sPlanck", "planck", "planck_wn"]
