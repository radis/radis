.. _label_bibliography:


Bibliography
------------

List of bibliographic references used in this project:

.. [CANTERA] D. G. Goodwin, H. K. Moffat, R. L. Speth, "Cantera: An Object-oriented Software
             Toolkit for Chemical Kinetics"", Thermodynamics, and Transport Processes,
             http://www.cantera.org, `doi:10.5281/zenodo.170284 <https://zenodo.org/record/170284#.XRIOno-xVEY>`__, 2017.

.. [RADIS-2018] E. Pannier, C. O. Laux, "RADIS: A Nonequilibrium Line-by-Line Radiative Code for CO2 and
                HITRAN-like database species", Journal of Quantitative Spectroscopy and Radiative Transfer
                (2018) `doi:10.1016/j.jqsrt.2018.09.027 <https://www.sciencedirect.com/science/article/pii/S0022407318305867?via%3Dihub>`__

.. [Spectral-Synthesis-Algorithm] D.C.M. van den Bekerom, E. Pannier,
                "A Discrete Integral Transform for Rapid Spectral Synthesis",
                Journal of Quantitative Spectroscopy and Radiative Transfer (2021)
                `doi:10.1016/j.jqsrt.2020.107476 <https://www.sciencedirect.com/science/article/abs/pii/S0022407320310049>`__

.. [Rothman-1998] L.S. Rothman, C.P. Rinsland, A. Goldman, S.T. Massie D.P. Edwards, J-M. Flaud,
                 A. Perrin, C. Camy-Peyret, V. Dana, J.-Y. Mandin, J. Schroeder, A. McCann,
                 R.R. Gamache, R.B. Wattson, K. Yoshino, K.V. Chance, K.W. Jucks, L.R. Brown,
                 V. Nemtchinov, P. Varanasi "The Hitran Molecular Spectroscopic Database
                 and Hawks (Hitran Atmospheric Workstation): 1996 Edition",
                 Journal of Quantitative Spectroscopy and Radiative Transfer 60 (1998)
                 665 - 710, `doi:10.1016/S0022-4073(98)00078-8 <https://www.sciencedirect.com/science/article/abs/pii/S0022407398000788?via%3Dihub>`__

.. [TIPS-2020] R.R Gamache, B. Vispoel, M. Rey, A. Nikitin, V. Tyuterev, O. Egorov, I.E Gordon, V. Boudon,
                "Total internal partition sums for the HITRAN2020 database",
                Journal of Quantitative Spectroscopy and Radiative Transfer 271 (2021)
                `doi:10.1016/j.jqsrt.2021.107713 <https://www.sciencedirect.com/science/article/abs/pii/S0022407321002065?via%3Dihub>`__

Line Databases
--------------

Reference of supported line databases:

.. [HITRAN-2016] I. Gordon, L. Rothman, C. Hill, R. Kochanov, Y. Tan, P. Bernath, V. Boudon, A. Campargue,
                 B. Drouin, J. M. Flaud, R. Gamache, J. Hodges, V. Perevalov, K. Shine, M.-a. Smith,
                 The HITRAN2016 Molecular Spectroscopic Database, Journal of Quantitative Spectroscopy and Radiative
                 Transfer 6 (38) (2017) 3–69, ISSN 00224073, `doi:10.1016/j.jqsrt.2017.06.038 <https://www.sciencedirect.com/science/article/pii/S0022407317301073>`__.

The latest HITRAN database version is automatically downloaded if using ``databank='hitran'``.

.. [HITEMP-2010] L. S. Rothman, I. E. Gordon, R. J. Barber, H. Dothe, R. R. Gamache, A. Goldman, V. I. Perevalov,
                 S. A. Tashkun, J. Tennyson, HITEMP, the high-temperature molecular spectroscopic database,
                 Journal of Quantitative Spectroscopy and Radiative Transfer 111 (15) (2010)
                 2139–2150, ISSN 00224073, `doi:10.1016/j.jqsrt.2010.05.001 <https://www.sciencedirect.com/science/article/pii/S002240731000169X>`__.

The latest HITEMP database version is automatically downloaded if using ``databank='hitran'``.

.. [CDSD-4000] S. A. Tashkun, V. I. Perevalov, CDSD-4000: High-resolution, high-temperature carbon dioxide
               spectroscopic databank, Journal of Quantitative Spectroscopy and Radiative Transfer 112 (9) (2011)
               1403–1410, ISSN 00224073, `doi:10.1016/j.jqsrt.2011.03.005 <https://www.sciencedirect.com/science/article/pii/S0022407311001154>`__

.. [ExoMol-2020] Tennyson et al., The 2020 release of the ExoMol database: Molecular line lists for
                exoplanet and other hot atmospheres, Journal of Quantitative Spectroscopy and Radiative Transfer 255,
                (2020), 107228,  `doi:10.1016/j.jqsrt.2020.107228 <https://www.sciencedirect.com/science/article/abs/pii/S002240732030491X>`__

For download and configuration of line databases, see the :ref:`Line Databases section <label_line_databases>`


Tools Used Within RADIS
-----------------------

For data retrieval :

.. [HAPI] `HAPI: The HITRAN Application Programming Interface <http://hitran.org/hapi>`_
          R. Kochanov, I. Gordon, L. Rothman, P. Wcisło, C. Hill, J. Wilzewski,
          "HITRAN Application Programming Interface (HAPI):
          A comprehensive approach to working with spectroscopic data", Journal of Quantitative Spectroscopy
          and Radiative Transfer 177 (2016) 15–30, ISSN 00224073, `doi:10.1016/j.jqsrt.2016.03.005 <https://www.researchgate.net/publication/297682202_HITRAN_Application_Programming_Interface_HAPI_A_comprehensive_approach_to_working_with_spectroscopic_data>`__.

.. [Astroquery] astroquery: An Astronomical Web-querying Package in Python 10.3847/1538-3881/aafc33
