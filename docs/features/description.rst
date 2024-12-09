Introduction
===========

Written as a general purpose radiative solver, the code is built around the spectroscopy databases
(`HITRAN/HITEMP <https://hitran.org>`_, `GEISA <https://geisa.aeris-data.fr/>`_, `ExoMol <https://www.exomol.com/>`_, etc.)
for molecules in their electronic ground state. Energy
levels are read from tabulated databases or calculated from Dunham developments.
Boltzmann, Treanor, and state specific vibrational distributions can be generated.
Spectra at thermal equilibrium can be computed for all species (:py:data:`~radis.db.MOLECULES_LIST_EQUILIBRIUM`).
and non-equilibrium spectra can be computed for |CO2| and CO
(:py:data:`~radis.db.MOLECULES_LIST_NONEQUILIBRIUM`).

To fit experimental spectra, RADIS includes a
:py:class:`~radis.tools.line_survey.LineSurvey` tool, an
interface with a look-up :py:class:`~radis.tools.database.SpecDatabase`
to improve fitting convergence times, and a :ref:`multi-slab module <label_los_index>`
with a radiative transfer equation solver to reproduce line-of-sight
experiments. :ref:`Validation cases <label_dev_test>` against existing
spectral codes and experimental results from various plasma sources are included [RADIS-2018]_.

.. |CO2| replace:: CO\ :sub:`2`

