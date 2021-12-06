# -*- coding: utf-8 -*-
"""
Doi  of some references used in this work.

See :ref:`the Bibliography page <label_bibliography>` for more information
"""

doi = {
    # Calculation / algorithms
    "RADIS-2018": "10.1016/j.jqsrt.2018.09.027",
    "DIT-2020": "10.1016/j.jqsrt.2020.107476",
    # Line Databases
    "HITRAN-1996": "10.1016/S0022-4073(98)00078-8",
    "HITRAN-2016": "10.1016/j.jqsrt.2017.06.038",
    "HITRAN-2020": "10.1016/j.jqsrt.2021.107949",
    "HITEMP-2010": "10.1016/j.jqsrt.2010.05.001",
    "CDSD-4000": "10.1016/j.jqsrt.2011.03.005",
    "ExoMol-2020": "10.1016/j.jqsrt.2020.107228",
    # Tabulated Partition functions
    "TIPS-2020": "10.1016/j.jqsrt.2021.107713",
    # Retrieval tools:
    "HAPI": "10.1016/j.jqsrt.2016.03.005",
    "Astroquery": "10.3847/1538-3881/aafc33",
    # Other
    "CANTERA": "10.5281/zenodo.170284",
    # Spectroscopic constants
    "Klarenaar-2017": "10.1088/1361-6595/aa902e",
    "Guelachvili-1983": "10.1016/0022-2852(83)90203-5",
}
"""dict: list of references used in RADIS calculations.

Print the BibTex entries used to compute your spectrum with :py:meth:`~radis.spectrum.spectrum.Spectrum.cite` ::

    s = calc_spectrum(...)
    s.cite()

"""
