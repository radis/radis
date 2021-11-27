Description
===========

Written as a general purpose radiative solver, the code is built around the [HITRAN-2020]_,
[HITEMP-2010]_ and [CDSD-4000]_ databases for molecules in their electronic ground state. Energy
levels are read from tabulated databases or calculated from Dunham developments.
Boltzmann, Treanor, and state specific vibrational distributions can be generated.
Thus far, |CO2|, CO are featured for non-equilibrium calculations
(:py:data:`~radis.db.MOLECULES_LIST_NONEQUILIBRIUM`),
and all species present in the HITRAN database are featured for equilibrium
calculations (:py:data:`~radis.db.MOLECULES_LIST_EQUILIBRIUM`).

To fit experimental spectra, RADIS includes a
:py:class:`~radis.tools.line_survey.LineSurvey` tool, an
interface with a look-up :py:class:`~radis.tools.database.SpecDatabase`
to improve fitting convergence times, and a :ref:`multi-slab module <label_los_index>`
with a radiative transfer equation solver to reproduce line-of-sight
experiments. :ref:`Validation cases <label_dev_test>` against existing
spectral codes and experimental results from various plasma sources are included [RADIS-2018]_.

.. |CO2| replace:: CO\ :sub:`2`

