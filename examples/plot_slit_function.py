# -*- coding: utf-8 -*-
"""
=============
Slit Function
=============

Quickly load and plot a slit function with :py:func:`~radis.tools.slit.plot_slit`

"""


from radis import plot_slit
from radis.test.utils import getTestFile

my_slit = getTestFile("slitfunction.txt")  # for the example here
plot_slit(my_slit, wunit="nm")
