# -*- coding: utf-8 -*-
"""
.. _example_hitran_full_range:

===============================
Calculate a full range spectrum
===============================

We compute the transmittance of an homogeneous, 1-km atmosphere layer
at 300 K and 1 atm, with 420 ppm CO2 and 2% H2O.

This example makes use of a lot of the RADIS optimization and convenience functions :

- The database is retrieved from the latest HITRAN data using [HAPI]_, parsed
  with RADIS and stored to memory-mapping formats. After the first calls,
  retrieving the database will be very fast

- The RADIS algorithm includes a sparse wavenumber implementation, which is
  activated by default in the following example where the waverange is very
  large. This behavior can be changed by setting the `SPARSE_WAVENUMBER` key
  of the :py:attr:`radis.config` dictionary, or of the ~/radis.json user file.

- The wavenumber grid is automatically computed using the `wstep='auto'` parameter.

- A linestrength cutoff reduces the initial number of lines.

Eventually, the full-range spectrum (about 4.8M points) is computed within a minute,
with about 140k lines resolved, i.e. a performance of about 1e10 gridpoints*lines/s.

Extra optimization could be achieved by cutting the spectrum in small intervals,
but requires an a-priori knowledge of the absorption & emission ranges of the
spectrum. See the :ref:`Calculate a large spectrum by part <example_large_range_by_part>`


"""


#%%

import astropy.units as u

from radis import calc_spectrum

s = calc_spectrum(
    wmin=0.5 * u.um,
    wmax=15 * u.um,  # cm-1
    mole_fraction={"CO2": 420e-6, "H2O": 0.02},
    isotope="1,2,3",
    pressure=1.01325,  # bar
    Tgas=300,  # K
    path_length=1e5,  # 1 km in cm
    verbose=2,
    databank="hitran",
    wstep="auto",
)


#%%
# Plot low and high resolution spectra on the same graph :
#
s.apply_slit(10, "nm")
import matplotlib.pyplot as plt

plt.figure(figsize=(16, 6))
s.plot("transmittance_noslit", wunit="nm", color="k", alpha=0.1, nfig="same")
s.plot("transmittance", wunit="nm", color="k", nfig="same")


#%%
# Print some details about the computed spectrum :
# We could also simply use `print(s)`

print("Number of grid points: ", len(s))
print(
    "Number of lines: ", s.c["total_lines"]
)  # some are discarded by linestrength cutoff
print("Number of lines: ", s.c["lines_calculated"])
print("Whether the sparse wavenumber algorithm was activated: ", s.c["sparse_ldm"])
print(f"Lineshape truncation used: {s.c['truncation']:.1f}cm-1")
print(f"Total calculation time: {s.c['calculation_time']:.1f}s")
print(
    f"Gridpoints * lines/s : {len(s)*s.c['lines_calculated']/s.c['calculation_time']:.1e}"
)
