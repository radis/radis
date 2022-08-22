# -*- coding: utf-8 -*-
"""
==================================
See populations of computed levels
==================================

Spectrum methods gives you access to the populations of levels used in a
nonequilibrium spectrum :

- :py:meth:`~radis.spectrum.spectrum.Spectrum.get_populations` provides the populations
of all levels of the molecule, used to compute the nonequilibrium partition function

- :py:attr:`~radis.spectrum.spectrum.Spectrum.lines` returns all informations
of the visible lines, in particular the populations of the upper and lower
levels of absorbing/emitting lines in the spectrum.

"""

from astropy import units as u

#%%
# We first compute a noneq CO spectrum.
# Notice the keys `export_lines=True` and `export_populations=True`
from radis import calc_spectrum

s = calc_spectrum(
    1900 / u.cm,
    2300 / u.cm,
    molecule="CO",
    isotope="1",
    pressure=1.01325 * u.bar,
    Tvib=3000 * u.K,
    Trot=1000 * u.K,
    mole_fraction=0.1,
    path_length=1 * u.cm,
    databank="hitran",  # or 'hitemp', 'geisa', 'exomol'
    export_lines=True,
    export_populations="rovib",
)

#%%
# Below we plot the populations
# We could also have used :py:meth:`~radis.spectrum.spectrum.Spectrum.plot_populations` directly

pops = s.get_populations("CO")["rovib"]

#%%
# Plot populations :
# Here we plot population fraction vs rotational number of the 2nd vibrational level ``v==2``,
# for all levels as well as for levels of lines visible in the spectrum

import matplotlib.pyplot as plt

pops.query("v==2").plot(
    "j",
    "n",
    label="populations used to compute\n partition functions",
    style=".-",
    markersize=2,
)
s.lines.query("vl==2").plot(
    "jl",
    "nl",
    ax=plt.gca(),
    label="populations of visible\n absorbing lines",
    kind="scatter",
    color="r",
)
plt.xlim((0, 80))
plt.legend()


#%%
# Print all levels information
print(f"{len(pops)} levels ")
print(pops)
print(pops.columns)


#%%
# Print visible lines information :
print(f"{len(s.lines)} visible lines in the spectrum")
print(s.lines)
print(s.lines.columns)

# %%
# We can also look at unique levels; for instance sorting by quantum numbers ``v, J``
# of the lower level : ``vl, jl``

print(
    len(s.lines.groupby(["vl", "jl"])),
    " visible absorbing levels (i.e. unique 'vl', 'jl' ",
)


#%%
#
# :py:meth:`~radis.spectrum.spectrum.Spectrum.line_survey` is also a convenient way to explore populatinos and other line
# parameters, for instance ::
#
#     s.line_survey(overlay="absorbance", barwidth=0.001, lineinfo="all")
#
#
# .. image:: https://user-images.githubusercontent.com/16088743/185366359-c6ecff54-95f2-4016-b762-649a7b4a9a1e.png
#     :alt: line survey output
#
