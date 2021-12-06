.. _label_dev_architecture:

Architecture
============

The RADIS modules are organized with the following flow chart

.. image:: RADIS_flow_chart.*
    :alt: https://radis.readthedocs.io/en/latest/_images/RADIS_flow_chart.svg
    :scale: 100 %

-------------------------------------------------------------------------

The upper part shows the successive calculation steps of the Line-by-Line module.
These steps appear clearly in the source code of the
:py:meth:`~radis.lbl.factory.SpectrumFactory.eq_spectrum` and
:py:meth:`~radis.lbl.factory.SpectrumFactory.non_eq_spectrum` methods of the
:py:class:`~radis.lbl.factory.SpectrumFactory`.
See details at the end of this file.

Methods are written in Factory objects inherited with the following scheme:

:py:class:`~radis.lbl.loader.DatabankLoader` > :py:class:`~radis.lbl.base.BaseFactory` >
:py:class:`~radis.lbl.broadening.BroadenFactory` > :py:class:`~radis.lbl.bands.BandFactory` >
:py:class:`~radis.lbl.factory.SpectrumFactory`


-------------------------------------------------------------------------

The Input Conditions in the left part, and the Computation Parameters on the right part,
are the input parameters of the different RADIS front-ends:

- :py:func:`~radis.lbl.calc.calc_spectrum` for the simple cases.
- :py:class:`~radis.lbl.factory.SpectrumFactory` with :py:meth:`~radis.lbl.factory.SpectrumFactory.eq_spectrum`
  and :py:meth:`~radis.lbl.factory.SpectrumFactory.non_eq_spectrum` for the other cases.
  GPU calculation can be done with :py:meth:`~radis.lbl.factory.SpectrumFactory.eq_spectrum_gpu`


-------------------------------------------------------------------------

The Input databases are either automatically downloaded from [HITRAN-2020]_, or defined
locally in a :ref:`Configuration file <label_lbl_config_file>`

-------------------------------------------------------------------------


The bottom part includes the post-processing modules of RADIS, in particular:

- The various methods associated with the :py:class:`~radis.spectrum.spectrum.Spectrum` class.

- The :ref:`Line-of-Sight module <label_los_index>` module

- The :py:class:`~radis.tools.line_survey.LineSurvey` tool.

- The :py:class:`~radis.tools.database.SpecDatabase` tool.



-------------------------------------------------------------------------

Methods from the Flow Chart: this methods are called successively from the
:py:meth:`radis.lbl.factory.SpectrumFactory.eq_spectrum` and
:py:meth:`radis.lbl.factory.SpectrumFactory.non_eq_spectrum` methods.

- Line Database: methods of :py:class:`~radis.lbl.loader.DatabankLoader` :

    - :py:meth:`radis.lbl.loader.DatabankLoader.load_databank`
    - :py:meth:`radis.lbl.loader.DatabankLoader.init_databank`
    - :py:meth:`radis.lbl.loader.DatabankLoader.fetch_databank`

- Partition functions: methods of :py:class:`~radis.levels.partfunc.RovibParFuncTabulator`
  and :py:class:`~radis.levels.partfunc.RovibParFuncCalculator` :

    - :py:meth:`radis.levels.partfunc.RovibParFuncTabulator.at`
    - :py:meth:`radis.levels.partfunc.RovibParFuncCalculator.at`
    - :py:meth:`radis.levels.partfunc.RovibParFuncCalculator.at_noneq`
    - :py:meth:`radis.levels.partfunc.RovibParFuncCalculator.at_noneq_3Tvib`

- Populations: methods of :py:class:`~radis.lbl.base.BaseFactory` :

    - :py:meth:`radis.lbl.base.BaseFactory.calc_populations_eq`
    - :py:meth:`radis.lbl.base.BaseFactory.calc_populations_noneq`

- Line Intensities: methods of :py:class:`~radis.lbl.base.BaseFactory` :

    - :py:meth:`radis.lbl.base.BaseFactory.calc_linestrength_eq`
    - :py:meth:`radis.lbl.base.BaseFactory._calc_linestrength_noneq`
    - :py:meth:`radis.lbl.base.BaseFactory._calc_emission_integral`

- Line Positions:  methods of :py:class:`~radis.lbl.base.BaseFactory` :

    - :py:meth:`radis.lbl.base.BaseFactory.calc_lineshift`

- Reduced line set: methods of :py:class:`~radis.lbl.base.BaseFactory` :

    - :py:meth:`radis.lbl.base.BaseFactory._cutoff_linestrength`

- Voigt Broadening: methods of :py:class:`~radis.lbl.broadening.BroadenFactory` :

    - :py:func:`radis.lbl.broadening.voigt_broadening_FWHM`
    - :py:func:`radis.lbl.broadening.voigt_lineshape`
    - :py:func:`radis.lbl.broadening._whiting`
    - :py:func:`radis.lbl.broadening._whiting_jit`
    - :py:meth:`radis.lbl.broadening.BroadenFactory._calc_broadening_FWHM`
    - :py:meth:`radis.lbl.broadening.BroadenFactory._add_voigt_broadening_FWHM`

- Pseudo-continuum: methods of :py:class:`~radis.lbl.broadening.BroadenFactory` :

    - :py:meth:`radis.lbl.broadening.BroadenFactory._find_weak_lines`
    - :py:meth:`radis.lbl.broadening.BroadenFactory._calculate_pseudo_continuum`
    - :py:meth:`radis.lbl.broadening.BroadenFactory._add_pseudo_continuum`

- Spectral densities k, j: methods of :py:class:`~radis.lbl.factory.SpectrumFactory` :

    - :py:meth:`radis.lbl.factory.SpectrumFactory.eq_spectrum`
    - :py:meth:`radis.lbl.factory.SpectrumFactory.non_eq_spectrum`

- RTE (1 slab): methods of :py:class:`~radis.lbl.factory.SpectrumFactory` :

    - :py:meth:`radis.lbl.factory.SpectrumFactory.eq_spectrum`
    - :py:meth:`radis.lbl.factory.SpectrumFactory.non_eq_spectrum`
