# -*- coding: utf-8 -*-
"""Print functions, with colors."""

import sys
from io import StringIO

import pandas as pd
from termcolor import colored

# %% Colored Print functions


def printm(*args, **kwargs):
    """Print in magenta.

    Typically use in RADIS for tests
    Could be used for Accuracy/Performance warnings."""

    string = _capture_print(*args, **kwargs)

    return print(colored(string, color="magenta"))


def printg(*args, **kwargs):
    """Print in green.

    Typically use in RADIS for verbose>=2 messages"""

    string = _capture_print(*args, **kwargs)

    return print(colored(string, color="green"))


def printr(*args, **kwargs):
    """Print in red.

    Typically use in RADIS to print errors that are corrected on-the-fly
    by the code"""

    string = _capture_print(*args, **kwargs)

    return print(colored(string, color="red"))


def _capture_print(*args, **kwargs):
    """Emulate print option but get output in a ```str``` variable instead of
    stdout."""

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
    """Print full Pandas series.

    From https://stackoverflow.com/questions/19124601/pretty-print-an-entire-pandas-series-dataframe
    """
    pd.set_option("display.max_rows", len(x))
    print(x)
    pd.reset_option("display.max_rows")


def get_print_full(x):
    """Same as print_full, but returns string."""

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


# %% UTF8 Upper script/ lower script
# https://stackoverflow.com/a/58612677/5622825

superscript_map = {
    "0": "⁰",
    "1": "¹",
    "2": "²",
    "3": "³",
    "4": "⁴",
    "5": "⁵",
    "6": "⁶",
    "7": "⁷",
    "8": "⁸",
    "9": "⁹",
    "a": "ᵃ",
    "b": "ᵇ",
    "c": "ᶜ",
    "d": "ᵈ",
    "e": "ᵉ",
    "f": "ᶠ",
    "g": "ᵍ",
    "h": "ʰ",
    "i": "ᶦ",
    "j": "ʲ",
    "k": "ᵏ",
    "l": "ˡ",
    "m": "ᵐ",
    "n": "ⁿ",
    "o": "ᵒ",
    "p": "ᵖ",
    "q": "۹",
    "r": "ʳ",
    "s": "ˢ",
    "t": "ᵗ",
    "u": "ᵘ",
    "v": "ᵛ",
    "w": "ʷ",
    "x": "ˣ",
    "y": "ʸ",
    "z": "ᶻ",
    "A": "ᴬ",
    "B": "ᴮ",
    "C": "ᶜ",
    "D": "ᴰ",
    "E": "ᴱ",
    "F": "ᶠ",
    "G": "ᴳ",
    "H": "ᴴ",
    "I": "ᴵ",
    "J": "ᴶ",
    "K": "ᴷ",
    "L": "ᴸ",
    "M": "ᴹ",
    "N": "ᴺ",
    "O": "ᴼ",
    "P": "ᴾ",
    "Q": "Q",
    "R": "ᴿ",
    "S": "ˢ",
    "T": "ᵀ",
    "U": "ᵁ",
    "V": "ⱽ",
    "W": "ᵂ",
    "X": "ˣ",
    "Y": "ʸ",
    "Z": "ᶻ",
    "+": "⁺",
    "-": "⁻",
    "=": "⁼",
    "(": "⁽",
    ")": "⁾",
}


subscript_map = {
    "0": "₀",
    "1": "₁",
    "2": "₂",
    "3": "₃",
    "4": "₄",
    "5": "₅",
    "6": "₆",
    "7": "₇",
    "8": "₈",
    "9": "₉",
    "a": "ₐ",
    "b": "♭",
    "c": "꜀",
    "d": "ᑯ",
    "e": "ₑ",
    "f": "բ",
    "g": "₉",
    "h": "ₕ",
    "i": "ᵢ",
    "j": "ⱼ",
    "k": "ₖ",
    "l": "ₗ",
    "m": "ₘ",
    "n": "ₙ",
    "o": "ₒ",
    "p": "ₚ",
    "q": "૧",
    "r": "ᵣ",
    "s": "ₛ",
    "t": "ₜ",
    "u": "ᵤ",
    "v": "ᵥ",
    "w": "w",
    "x": "ₓ",
    "y": "ᵧ",
    "z": "₂",
    "A": "ₐ",
    "B": "₈",
    "C": "C",
    "D": "D",
    "E": "ₑ",
    "F": "բ",
    "G": "G",
    "H": "ₕ",
    "I": "ᵢ",
    "J": "ⱼ",
    "K": "ₖ",
    "L": "ₗ",
    "M": "ₘ",
    "N": "ₙ",
    "O": "ₒ",
    "P": "ₚ",
    "Q": "Q",
    "R": "ᵣ",
    "S": "ₛ",
    "T": "ₜ",
    "U": "ᵤ",
    "V": "ᵥ",
    "W": "w",
    "X": "ₓ",
    "Y": "ᵧ",
    "Z": "Z",
    "+": "₊",
    "-": "₋",
    "=": "₌",
    "(": "₍",
    ")": "₎",
}
