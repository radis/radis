
**************************
Line-of-sight (LOS) module
**************************

This module takes several Spectrum objects as an input and combines then along the 
line-of-sight (SerialSlabs) or at the same spatial position (MergeSlabs), to 
reproduce line-of-sight experiments 


How to combine slabs?
=====================


Along the line-of-sight
-----------------------

Use the :func:`~radis.los.slabs.SerialSlabs` function: 

        >>> s1 = calc_spectrum(...)
        >>> s2 = calc_spectrum(...)
        >>> s3 = SerialSlabs(s1, s2)
        
        
At the same spatial position
----------------------------

Use the :func:`~radis.los.slabs.MergeSlabs` function:

Merge two spectra calculated with different species (true only if broadening
coefficient dont change much):

    >>> from radis import calc_spectrum, MergeSlabs
    >>> s1 = calc_spectrum(...)
    >>> s2 = calc_spectrum(...)
    >>> s3 = MergeSlabs(s1, s2)
    
Load a spectrum precalculated on several partial spectral ranges, for a same 
molecule (i.e, partial spectra are optically thin on the rest of the spectral 
range)

    >>> from radis import load_spec, MergeSlabs
    >>> spectra = []
    >>> for f in ['spec1.spec', 'spec2.spec', ...]:
    >>>     spectra.append(load_spec(f))
    >>> s = MergeSlabs(*spectra, resample_wavespace='full', out_of_bounds='transparent')
    >>> s.update()   # Generate missing spectral quantities
    >>> s.plot()
    
    
  