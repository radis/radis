# -*- coding: utf-8 -*-
"""
=============================
Load an experimental spectrum
=============================

Load an experimental spectrum stored as ``.spec`` (with units and metadata)
using :py:func:`~radis.tools.database.load_spec` (we could also have used
:py:func:`~radis.tools.database.plot_spec` directly !)

See more loading and post-processing functions on the :ref:`Spectrum page <label_spectrum>`.

"""


from radis.test.utils import getTestFile
from radis.tools.database import load_spec

my_file = getTestFile("CO2_measured_spectrum_4-5um.spec")  # for the example here

s = load_spec(my_file)
s.plot()

# Print all metadata:
print(s)

# Or retrieve an information :
print(s.crop(4160, 4200, "nm").max())
