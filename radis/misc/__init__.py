# -*- coding: utf-8 -*-
"""
Created on Tue May 26 11:52:15 2015

Erwan Pannier
EM2C, CentraleSupélec, 2015
CNRS UPR 288

"""

from __future__ import absolute_import

from .arrays import (norm, norm_on, scale_to, shift_array, calc_diff, 
                    find_nearest, find_first, autoturn, centered_diff, 
                    evenly_distributed, bining, count_nans, array_allclose,
                    nantrapz, logspace)
from .basics import (key_max_val, exec_file, remove_duplicates, partition,
                    is_float, list_if_float, compare_dict, compare_lists, compare_paths,
                    merge_lists)
from .curve import (curve_add, curve_substract, curve_divide, curve_multiply,
                    curve_distance)
from .config import (getDatabankEntries,
                    getDatabankList)
from .debug import export
from .signal import resample, resample_even
from .progress_bar import ProgressBar
from .utils import (getProjectRoot, NotInstalled, DatabankNotFound)