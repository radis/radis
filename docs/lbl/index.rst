
*************************
Line-by-line (LBL) module
*************************

This is the core of RADIS: it calculates the spectral densities for a homogeneous
slab of gas, and returns a :py:class:`~radis.spectrum.spectrum.Spectrum` object. 

The Spectrum Factory
--------------------

The :class:`~radis.lbl.factory.SpectrumFactory` allows you to 
load your own line database with :meth:`~radis.lbl.loader.DatabankLoader.load_databank`, 
and then calculate several spectra in batch using 
:meth:`~radis.lbl.factory.SpectrumFactory.eq_spectrum` and 
:meth:`~radis.lbl.factory.SpectrumFactory.non_eq_spectrum` ::

    from radis import SpectrumFactory
    sf = SpectrumFactory(wavelength_min=4150, 
                         wavelength_max=4400,
                         path_length=1,             # cm
                         pressure=0.020,            # bar
                         molecule='CO2',
                         isotope='1,2', 
                         cutoff=1e-25,              # cm/molecule  
                         broadening_max_width=10,   # cm-1
                         )
    sf.load_databank('CDSD-HITEMP')        # this database must be defined in ~/.radis
    s1 = sf.eq_spectrum(Tgas=300)
    s2 = sf.eq_spectrum(Tgas=2000)
    s3 = sf.non_eq_spectrum(Tvib=2000, Trot=300)





