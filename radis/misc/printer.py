# -*- coding: utf-8 -*-
"""
Print functions, with colors
"""

from __future__ import print_function
from termcolor import colored
from six import StringIO
import sys

# %% Colored Print functions


def printm(*args, **kwargs):
    ''' Print in magenta'''

    string = _capture_print(*args, **kwargs)

    return print(colored(string, color='magenta'))


def printg(*args, **kwargs):
    ''' Print in green'''

    string = _capture_print(*args, **kwargs)

    return print(colored(string, color='green'))


def printr(*args, **kwargs):
    ''' Print in red'''

    string = _capture_print(*args, **kwargs)

    return print(colored(string, color='red'))


def _capture_print(*args, **kwargs):
    ''' Emulate print option but get output in a ```str``` variable instead of stdout '''

    # Change the output to capture the string instead of sending it to the console
    old_stdout = sys.stdout
    sys.stdout = newstdout = StringIO()     # capture all print
    try:
        print(*args, **kwargs)
        string = newstdout.getvalue()              # Get string output
    except:
        raise
    finally:
        sys.stdout = old_stdout

    # discard line return character (will be added again when we print ``string`` )
    return string[:-1]
