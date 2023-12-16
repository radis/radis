.. _label_examples:

========
Examples
========

.. toctree::
   :maxdepth: 3

   examples

Many other examples scripts are available on the `radis-examples project <https://github.com/radis/radis-examples>`__.


CH4 Full Spectrum Benchmark
===========================

Here we reproduce the full spectrum (0.001 - 11500 cm-1) of Methane for a ``broadening_max_width``
corresponding to about 50 HWHMs, as in the Benchmark case of [HAPI]_, Table 7, Methane_III,
also featured in the [RADIS-2018]_ article ::

    from radis import SpectrumFactory

    benchmark_line_brd_ratio = 50    # “WavenumberWingHW”/HWHMs
    dnu = 0.01         # step in HAPI Benchmark article
    molecule = 'CH4'
    wavenum_min = 0.001
    wavenum_max = 11505
    pressure_bar = 1.01315
    T = 296
    isotopes = [1, 2, 3, 4]

    sf = SpectrumFactory(wavenum_min=wavenum_min,
                         wavenum_max=wavenum_max,
                         isotope=isotopes,  #'all',
                         verbose=2,
                         wstep=dnu,     # depends on HAPI benchmark.
                         cutoff=1e-23,
                         broadening_max_width=5.73,  # Corresponds to WavenumberWingHW/HWHM=50 in HAPI
                         molecule=molecule,
                         optimization=None,
                         )
    sf.fetch_databank('astroquery')

    s = sf.eq_spectrum(Tgas=T, pressure=pressure_bar)
    s.plot()


The comparison in terms of performance with HAPI can be found in the ``radis/test/benchmark/radis_vs_hapi_CH4_full_spectrum.py``
case::

    cd radis
    python radis/test/benchmark/radis_vs_hapi_CH4_full_spectrum.py

Using the different :ref:`Performance <label_lbl_performance>` optimizations available in RADIS,
the calculation is typically 100 times faster in RADIS::

    >>> Calculated with HAPI in 157.41s
    >>> Calculated with RADIS in 1.65s


