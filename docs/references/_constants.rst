.. _label_db_spectroscopic_constants:

Spectroscopic constants
=======================

CO2
---

Version 0.9 uses uses the CO2 spectrosopic coefficients of [Klarenaar2017]_ ,
Table 2,3 and the references therein. ⚠️ These constants have been compiled for Treanor distributions
and may not be suited for Boltzmann distributions.

.. include:: ../../radis/db/CO2/molecules_data.json
    :literal:
    :end-line: 20

See the full list of constants of the public version
`as text <https://raw.githubusercontent.com/radis/radis/master/radis/db/CO2/molecules_data.json>`__
or `on GitHub <https://github.com/radis/radis/blob/master/radis/db/CO2/molecules_data.json>`__

In version 1.0, RADIS will use the spectroscopic constants of [Suzuki1968]_

CO
--

CO uses the Dunham coefficients of [Guelachvili1983]_

.. include:: ../../radis/db/CO/molecules_data.json
    :literal:
    :end-line: 20

See the full list of constants of the public version
`as text <https://raw.githubusercontent.com/radis/radis/master/radis/db/CO/molecules_data.json>`__
or `on GitHub <https://github.com/radis/radis/blob/master/radis/db/CO/molecules_data.json>`__



You can use :ref:`your own set of spectroscopic constants <label_lbl_custom_constants>`,
or precompute energy levels and use them directly (see the
`Energy level database <https://radis.readthedocs.io/en/latest/lbl/lbl.html#label-lbl-config-file>`__).


References
==========

.. [Klarenaar2017] B. L. M. Klarenaar, R. Engeln, D. C. M. van den Bekerom, M. C. M. van de Sanden,
                   A. S. MorilloCandas, O. Guaitella, "Time evolution of vibrational temperatures
                   in a CO2 glow discharge measured with infrared absorption spectroscopy",
                   *Plasma Sources Science and Technology* 26 (11) (2017) 115008, ISSN 1361-6595,
                   `doi:10.1088/1361-6595/aa902e <https://iopscience.iop.org/article/10.1088/1361-6595/aa902e>`__.

.. [Suzuki1968] I. Suzuki, "General anharmonic force constants of carbon dioxide"
                *Journal of Molecular Spectroscopy* 25 479-500 ISSN 00222852
                `doi:10.1016/S0022-2852(68)80018-9 <https://www.sciencedirect.com/science/article/pii/S0022285268800189>`__

.. [Guelachvili1983] G.Guelachvili, D.de Villeneuve R.Farrenq, W.Urban, J.Verges,
                    Dunham coefficients for seven isotopic species of CO
                    `doi:10.1016/0022-2852(83)90203-5 <https://www.sciencedirect.com/science/article/pii/0022285283902035>`__
