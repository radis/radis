# -*- coding: utf-8 -*-
"""
Print functions, with colors
"""

from __future__ import absolute_import, unicode_literals, print_function, division
from termcolor import colored
from six import StringIO
import sys
import pandas as pd

# %% Colored Print functions


def printm(*args, **kwargs):
    """ Print in magenta"""

    string = _capture_print(*args, **kwargs)

    return print(colored(string, color="magenta"))


def printg(*args, **kwargs):
    """ Print in green"""

    string = _capture_print(*args, **kwargs)

    return print(colored(string, color="green"))


def printr(*args, **kwargs):
    """ Print in red"""

    string = _capture_print(*args, **kwargs)

    return print(colored(string, color="red"))


def _capture_print(*args, **kwargs):
    """ Emulate print option but get output in a ```str``` variable instead of stdout """

    # Change the output to capture the string instead of sending it to the console
    old_stdout = sys.stdout
    sys.stdout = newstdout = StringIO()  # capture all print
    try:
        print(*args, **kwargs)
        string = newstdout.getvalue()  # Get string output
    except:
        raise
    finally:
        sys.stdout = old_stdout

    # discard line return character (will be added again when we print ``string`` )
    return string[:-1]


# %% Print in pandas


def print_full(x):
    """ Print full Pandas series. 
    From https://stackoverflow.com/questions/19124601/pretty-print-an-entire-pandas-series-dataframe
    """
    pd.set_option("display.max_rows", len(x))
    print(x)
    pd.reset_option("display.max_rows")


def get_print_full(x):
    """ Same as print_full, but returns string"""

    # Change the output to capture the string instead of sending it to the console
    old_stdout = sys.stdout
    sys.stdout = newstdout = StringIO()  # capture all print
    try:
        pd.set_option("display.max_rows", len(x))
        print(x)
        string = newstdout.getvalue()  # Get string output
    except:
        raise
    finally:
        sys.stdout = old_stdout
        pd.reset_option("display.max_rows")

    return string[:-1]
