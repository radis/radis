# -*- coding: utf-8 -*-
"""Parsers for various databases

-------------------------------------------------------------------------------

"""

from .exomol import fetch_exomol
from .geisa import fetch_geisa
from .hitemp import fetch_hitemp
from .chianti import load_chianti
from .hitran import fetch_hitran
from .query import fetch_astroquery
__all__ = [
    "fetch_hitemp",
    "fetch_hitran",
    "fetch_exomol",
    "fetch_astroquery",
    "fetch_geisa",
    "load_chianti",
]

_lazy_imports = {
    "fetch_exomol": ".exomol",
    "fetch_geisa": ".geisa",
    "fetch_hitemp": ".hitemp",
    "fetch_hitran": ".hitran",
    "fetch_astroquery": ".query",
}


def __getattr__(name):
    if name in _lazy_imports:
        import importlib

        module = importlib.import_module(_lazy_imports[name], __name__)
        attr = getattr(module, name)
        globals()[name] = attr  # cache for subsequent access
        return attr
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(set(list(globals().keys()) + list(__all__)))
