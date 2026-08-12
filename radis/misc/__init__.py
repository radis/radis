# -*- coding: utf-8 -*-
"""Misc. and support functions"""

import importlib

__all__ = [
    "array_allclose",
    "autoturn",
    "bining",
    "calc_diff",
    "centered_diff",
    "count_nans",
    "evenly_distributed",
    "find_first",
    "find_nearest",
    "is_sorted",
    "is_sorted_backward",
    "logspace",
    "nantrapz",
    "norm",
    "norm_on",
    "scale_to",
    "compare_dict",
    "compare_lists",
    "compare_paths",
    "exec_file",
    "is_float",
    "key_max_val",
    "list_if_float",
    "make_folders",
    "merge_lists",
    "partition",
    "remove_duplicates",
    "getDatabankEntries",
    "getDatabankList",
    "curve_add",
    "curve_distance",
    "curve_divide",
    "curve_multiply",
    "curve_substract",
    "DatabaseProgressPrinter",
    "get_progress_printer",
    "export",
    "ProgressBar",
    "resample",
    "resample_even",
    "DatabankNotFound",
    "NotInstalled",
    "getProjectRoot",
]

_EXPORT_TO_MODULE = {
    "array_allclose": "arrays",
    "autoturn": "arrays",
    "bining": "arrays",
    "calc_diff": "arrays",
    "centered_diff": "arrays",
    "count_nans": "arrays",
    "evenly_distributed": "arrays",
    "find_first": "arrays",
    "find_nearest": "arrays",
    "is_sorted": "arrays",
    "is_sorted_backward": "arrays",
    "logspace": "arrays",
    "nantrapz": "arrays",
    "norm": "arrays",
    "norm_on": "arrays",
    "scale_to": "arrays",
    "compare_dict": "basics",
    "compare_lists": "basics",
    "compare_paths": "basics",
    "exec_file": "basics",
    "is_float": "basics",
    "key_max_val": "basics",
    "list_if_float": "basics",
    "make_folders": "basics",
    "merge_lists": "basics",
    "partition": "basics",
    "remove_duplicates": "basics",
    "getDatabankEntries": "config",
    "getDatabankList": "config",
    "curve_add": "curve",
    "curve_distance": "curve",
    "curve_divide": "curve",
    "curve_multiply": "curve",
    "curve_substract": "curve",
    "DatabaseProgressPrinter": "database_progress",
    "get_progress_printer": "database_progress",
    "export": "debug",
    "ProgressBar": "progress_bar",
    "resample": "signal",
    "resample_even": "signal",
    "DatabankNotFound": "utils",
    "NotInstalled": "utils",
    "getProjectRoot": "utils",
}


def __getattr__(name):
    module_name = _EXPORT_TO_MODULE.get(name)
    if module_name is None:
        raise AttributeError(f"module '{__name__}' has no attribute '{name}'")

    module = importlib.import_module(f".{module_name}", __name__)
    value = getattr(module, name)
    globals()[name] = value
    return value


def __dir__():
    return sorted(set(globals()) | set(__all__))
