# -*- coding: utf-8 -*-
"""Misc. and support functions
"""


from .arrays import (
    array_allclose,
    autoturn,
    bining,
    calc_diff,
    centered_diff,
    count_nans,
    evenly_distributed,
    find_first,
    find_nearest,
    is_sorted,
    is_sorted_backward,
    logspace,
    nantrapz,
    norm,
    norm_on,
    scale_to,
)
from .basics import (
    compare_dict,
    compare_lists,
    compare_paths,
    exec_file,
    is_float,
    key_max_val,
    list_if_float,
    make_folders,
    merge_lists,
    partition,
    remove_duplicates,
)
from .config import getDatabankEntries, getDatabankList
from .curve import (
    curve_add,
    curve_distance,
    curve_divide,
    curve_multiply,
    curve_substract,
)
from .debug import export
from .progress_bar import ProgressBar
from .signal import resample, resample_even
from .utils import DatabankNotFound, NotInstalled, getProjectRoot
