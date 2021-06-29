# -*- coding: utf-8 -*-
"""
.. _example_cite:

========================
Cite all references used
========================

RADIS is built on the shoulders of many state-of-the-art packages and databases. If using RADIS
for your work, **cite all of them that made it possible**.

Starting from 0.9.30, you can retrieve the bibtex entries of all papers and
references that contribute to the calculation of a Spectrum, with the
:py:meth:`~radis.spectrum.spectrum.Spectrum.cite` method.

Example
-------

Below, we compute a non-equilibrium spectrum. The :py:meth:`~radis.spectrum.spectrum.Spectrum.cite`
method returns the computational algorithm used, the line database,
the spectroscopic constants used to compute rovibrational energies, and the data retrieval tools.

See Also
--------
:py:class:`~radis.tools.track_ref.RefTracker`

"""


from radis import calc_spectrum

s = calc_spectrum(
    1900,
    2300,  # cm-1
    molecule="CO",
    isotope="1,2,3",
    pressure=1.01325,  # bar
    Tvib=2000,  #
    Trot=300,
    mole_fraction=0.1,
    path_length=1,  # cm
    databank="hitran",
)

s.cite()
