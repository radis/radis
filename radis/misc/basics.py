# -*- coding: utf-8 -*-
"""Created on Wed Nov  5 12:59:37 2014.

@author: Erwan
Small functions used in other procedures
-------------------------------------------------------------------------------
"""


import os
import sys
from io import StringIO
from itertools import filterfalse, tee
from os.path import abspath, join, normcase, normpath

import numpy as np
import pandas as pd

verbose = True

# %%
# ==============================================================================
#  OS Functions
# ==============================================================================


def make_folders(path, folders):
    """Make folders if not there
    Parameters
    ----------
    path: str
        where to create folders
    folders: list or str
        folders to create
    """
    if type(folders) is str:
        folders = [folders]
    for folder in folders:
        os.makedirs(
            join(path, folder), exist_ok=True
        )  # makedirs create folders recursively


# %%
# ==============================================================================
# Basic Functions
# ==============================================================================


def all_in(keys, L):
    """Returns whether all items in keys are in list L."""
    return all([k in L for k in keys])


def in_all(key, list_of_list):
    """Returns true if key is in all lists."""
    return all([key in L for L in list_of_list])


def any_in(keys, L):
    """Returns whether any of the items in keys are in list L."""
    return any([k in L for k in keys])


def key_max_val(d):
    """Return the dictionary key with max value."""
    v = list(d.values())
    k = list(d.keys())
    return k[v.index(max(v))]


def exec_file(afile, globalz=None, localz=None):
    with open(afile, "r") as fh:
        exec(fh.read(), globalz, localz)


def remove_duplicates(l):
    """Remove duplicates from a list, without changing the order.

    Note that if the order doesn't matter you could just do set(l)
    """

    l1 = []
    for e in l:
        if not e in l1:
            l1.append(e)
    return l1


def partition(pred, iterable):
    """Use a predicate to partition entries into false entries and true entries
    Returns
    -------
    Returns two lists: positive, and negative
    Example
    -------
     >>> partition(is_odd, range(10))
     --> [0 2 4 6 8], [1 3 5 7 9]
    """
    t1, t2 = tee(iterable)
    return list(filter(pred, t2)), list(filterfalse(pred, t1))


# Just for performance comparison. Below: naive implementation of partition
# %timeit partition(lambda x:x>0.5, np.random.rand(10))
# ... partition         2.67 µs per loop
# ... partition_naive   28.5 µs per loop
# def partition_naive(pred, iterable):
#    trues = []
#    falses = []
#    for item in iterable:
#        if pred(item):
#            trues.append(item)
#        else:
#            falses.append(item)
#    return trues, falses


def str2bool(s):
    """Used to convert Pandas columns in str type to boolean
    (note that by default bool("False")==True !)
    """
    if s in [True, "true", "True", 1, "1", "1.0"]:
        return True
    elif s in [False, "false", "False", 0, "0", "0.0"]:
        return False
    else:
        raise ValueError


# %% Compare / merge tools


def intersect(a, b):
    """Returns intersection of two dictionaries on values."""
    c = {}
    for k in set(a.keys()) & set(b.keys()):
        c[k] = a[k] if (a[k] == b[k]) else "N/A"
    return c


def compare_dict(
    d1,
    d2,
    verbose="if_different",
    compare_as_paths=[],
    compare_as_close=[],
    return_string=False,
    df1_str="Left",
    df2_str="Right",
    ignore_keys=[],
):
    """Returns ratio of equal keys [0-1]
    If verbose, also print all keys and values on 2 columns
    Parameters
    ----------
    d1, d2: dict
        two dictionaries to compare

    Other Parameters
    ----------------
    compare_as_paths: list of keys
        compare the values corresponding to given keys as path (irrespective of
        forward / backward slashes, or case )
    compare_as_close: list of keys
        compare with ``np.isclose(a,b)`` rather than ``a==b``
    verbose: boolean, or ``'if_different'``
        ``'if_different'`` means results will be shown only if there is a difference.
    return_string: boolean
        if ``True``, returns message instead of just printing it (useful in error messages)
        Default ``False``
    ignore_keys: list
        do not compare these keys

    Returns
    -------
    out: float [0-1]
        ratio of matching keys
    if ``return_string``:
    out, string: float [0-1], str
        ratio of matching keys and comparison message
    """

    old_stdout = sys.stdout
    sys.stdout = newstdout = StringIO()  # capture all print

    try:

        print("{0:15}{1:17}\t{2}".format("Key", df1_str, df2_str))
        print("-" * 40)
        all_keys = set(list(d1.keys()) + list(d2.keys()))
        all_keys = [k for k in all_keys if k not in ignore_keys]
        s = 0  # counter of all matching keys
        for k in all_keys:
            if k in d1 and k in d2:  # key in both dicts. Let's compare values
                # Deal with path case
                if k in compare_as_paths:
                    if not compare_paths(d1[k], d2[k]):
                        print("{0:15}{1}\t{2}".format(k, d1[k], d2[k]))
                    else:
                        s += 1
                # Deal with is close case
                elif k in compare_as_close:
                    if not np.isclose(d1[k], d2[k]):
                        print("{0:15}{1}\t{2}".format(k, d1[k], d2[k]))
                    else:
                        s += 1
                # Other cases
                else:
                    if d1[k] != d2[k]:
                        print("{0:15}{1}\t{2}".format(k, d1[k], d2[k]))
                    else:
                        s += 1
            elif k in d1:
                print("{0:15}{1}\tN/A".format(k, d1[k]))
            else:
                print("{0:15}N/A\t{1}".format(k, d2[k]))
        print("-" * 40)

        if len(all_keys) == 0:
            out = 1
        else:
            out = s / len(all_keys)

        # Get string output
        string = newstdout.getvalue()
        sys.stdout = old_stdout  # reset normal print

        # Output
        if verbose == True or (verbose == "if_different" and out != 1):
            print(string)

        if return_string:
            out = out, string

        return out
    except:
        raise
    finally:
        sys.stdout = old_stdout


def compare_lists(
    l1,
    l2,
    verbose="if_different",
    return_string=False,
    l1_str="Left",
    l2_str="Right",
    print_index=False,
):
    """Compare 2 lists of elements that may not be of the same length, irrespective
    of order. Returns the ratio of elements [0-1] present in both lists. If verbose,
    prints the differences

    Parameters
    ----------
    l1, l2: list-like
    verbose: boolean, or 'if_different'
        'if_different' means results will be shown only if there is a difference.
        function is called twice

    Other Parameters
    ----------------
    verbose: boolean, or ``'if_different'``
        ``'if_different'`` means results will be shown only if there is a difference.
    return_string: boolean
        if ``True``, returns message instead of just printing it (useful in error messages)
        Default ``False``

    Returns
    -------
    out: float [0-1]
        ratio of matching keys
    """

    old_stdout = sys.stdout
    sys.stdout = newstdout = StringIO()  # capture all print

    try:
        tab = "        " if print_index else ""
        print(tab + "{0:20}\t\t{1}".format(l1_str, tab + l2_str))
        print(tab + "-" * (44 + len(tab)))
        all_keys = set(list(l1) + list(l2))
        s = 0  # counter of all matching keys
        for i, k in enumerate(all_keys):
            if k in l1 and k in l2:  # key in both lists
                s += 1
            elif k in l1:
                l1_index_str = f"|#{l1.index(k):3}|  " if print_index else ""
                print(
                    "{0:20}\t\t{1}".format(
                        l1_index_str + "{0} ({1})".format(k, type(k)), tab + "N/A"
                    )
                )
            else:
                l2_index_str = f"|#{l2.index(k):3}|  " if print_index else ""
                print(
                    "{0:20}\t\t{1} ({2})".format(tab + "N/A", l2_index_str + k, type(k))
                )
        print(tab + "-" * (44 + len(tab)))

        if len(all_keys) == 0:
            out = 1
        else:
            out = s / len(all_keys)

        # Get string output
        string = newstdout.getvalue()
        sys.stdout = old_stdout  # reset normal print

        # Output
        if verbose == True or (verbose == "if_different" and out != 1):
            print(string)

        if return_string:
            out = out, string

        return out
    except:
        raise
    finally:
        sys.stdout = old_stdout


def stdpath(p):
    """Convert path p in standard path (irrespective of slash / backslash, or
    case)"""

    return normpath(normcase(abspath(p)))


def compare_paths(p1, p2):
    """Compare 2 paths p1 and p2."""
    return stdpath(p1) == stdpath(p2)


def merge_lists(lists):
    """Merge a list of lists and return a list with unique elements."""
    return list(set(sum([l for l in lists], [])))


# ==============================================================================
# %% Pandas specific
# ==============================================================================


def merge_rename_columns(df, columns1, columns2, merged_names):
    """Merge all columns under easier names. Only keep the useful ones
    Returns a new dataframe
    Parameters
    ----------
    df: pandas Dataframe
    columns1: list
        list of columns names
    columns2: list
        list of columns names, whose index match columns 1
    merged_names: list
        new names
    Example
    -------
    df = merge_rename_columns(df1, ['lvl_u', 'ju', 'Eu', 'nu', 'gu', 'grotu'],
                                   ['lvl_l', 'jl', 'El', 'nl', 'gl', 'grotl'],
                                   ['lvl',   'j',  'E',  'n',  'g',  'grot']
                                   )
    """

    assert all_in(columns1, list(df.keys()))
    assert all_in(columns2, list(df.keys()))

    df1 = df.loc[:, columns1]
    df2 = df.loc[:, columns2]
    df1.rename(
        columns={columns1[i]: merged_names[i] for i in range(len(merged_names))},
        inplace=True,
    )
    df2.rename(
        columns={columns2[i]: merged_names[i] for i in range(len(merged_names))},
        inplace=True,
    )
    df = pd.concat((df1, df2), ignore_index=True)

    return df.drop_duplicates()


def print_series(a):
    """Print a pandas series `a` , explicitely showing all rows."""

    for i, k in enumerate(a.keys()):
        print(k, "\t", a.values[0][i])


def transfer_metadata(df1, df2, metadata):
    """Transfer metadata between a DataFrame df1 and df2.

    For some reason metadata are sometimes not copied when a DataFrame is
    sliced or copied, even if they explicitely appear in the ``df.attrs``
    attribute. See https://github.com/pandas-dev/pandas/issues/28283

    Here we copy them back. Attributes can be :
        - keys of the ``df1.attrs`` dictionary
        - simple attributes of ``df1``, i.e., ``df1.X``

    Parameters
    ----------
    df1: pandas DataFrame
        copy from df1
    df2: pandas DataFrame
        copy to df2
    """

    for k in metadata:
        if not hasattr(df2, k):
            assert k not in df1.columns  # point is not to copy columns!
            if k in df1.attrs:  # Keys of attribute dictionary
                df2.attrs[k] = df1.attrs[k]
            else:  # Direct attributes of df1
                setattr(df2, k, getattr(df1, k))

    # @dev: refactor : we're updating the database to properly store values
    # either in columns either in the .attrs dict, but so they can always
    # be accessed with df.X


def expand_metadata(df, metadata):
    """Turn metadata from a float to a column.

    For some reason metadata are sometimes not copied when a DataFrame is
    sliced or copied, even if they explicitely figure in the df.attrs
    attribute. Here we add them as column before such operations.

    Parameters
    ----------

    df: pandas DataFrame
        ...

    Returns
    -------

    None:
        df modified in place
    """

    for k in metadata:
        if __debug__ and k not in df.attrs:
            from radis.misc.debug import printdbg

            printdbg("WARNING. {0} not in metadata: {1}".format(k, df.attrs))
        df[k] = getattr(df, k)


# %%
# ==============================================================================
# Types
# ==============================================================================


def list_if_float(a):
    if type(a) is list:
        return a
    else:
        return [a]


def flatten(*args):
    """Flatten list of lists of floats."""
    out = []
    for a in args:
        if is_list(a):
            out += list(a)
        else:
            out += [a]
    return out


def is_list(a):
    """Returns True if a has list-like type: list, np.array, tuple, set,
    etc.)"""
    return type(a) in [list, np.ndarray, tuple, set]


def is_float(a):
    """Returns True if a has float-like type: float, np.float64, np.int64,
    etc.)"""
    return type(a) in [float, np.float64, np.int32, np.float32, int, np.int64]


def is_range(a):
    from radis.tools.plot_tools import ParamRange

    return isinstance(a, ParamRange)


def is_number(s):
    """Return True if ``s`` is a number.

    Works for strings, floats, int, and is compatible with Python 2/3
    """

    try:
        float(s)
    except (ValueError, TypeError):
        return False
    else:
        return True


def to_str(a):
    if isinstance(a, bytes):
        return a.decode("utf-8")
    else:
        return a


def round_off(n):
    # Getting rounded off value (atleast order 3)
    for i in range(0, 10):
        val = round(n, 3 + i)
        if val == 0:
            continue
        return val
    return 0
