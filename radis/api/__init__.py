# -*- coding: utf-8 -*-
"""Common API for RADIS & `Exojax <https://github.com/HajimeKawahara/exojax>`__

-------------------------------------------------------------------------------

"""

from .cdsdapi import cdsd2df
from .hitranapi import hit2df

__all__ = [
    "cdsd2df",
    "hit2df",
]
