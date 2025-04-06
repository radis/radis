# -*- coding: utf-8 -*-
"""Created on Thu Aug 13 14:51:22 2020.

@author: erwan
"""

import numpy as np
import pandas as pd


def npy2df(keywords, verbose=True):
    """Convert a dictionary of Numpy arrays storing spectroscopic information
    to a Pandas dataframe.

    Parameters
    ----------
    keywords: dict

        npy2df({'wav':'PATH/TO/v0.npy',
          'int':'PATH/TO/int.npy',
          'Pshft':'PATH/TO/int.npy',
          'log_2gs':'PATH/TO/log_2gs.npy'
          'Tdpair':'PATH/TO/Tdpair.npy',
          'El':'PATH/TO/Tdpair.npy'}

    Examples
    --------
    ::

        npy2df({'wav':'PATH/TO/v0.npy',
          'int':'PATH/TO/int.npy',
          'Pshft':'PATH/TO/int.npy',
          'log_2gs':'PATH/TO/log_2gs.npy'
          'Tdpair':'PATH/TO/Tdpair.npy',
          'El':'PATH/TO/Tdpair.npy'}

    See definitions for instance in :py:data:`~radis.api.hitranapi.column_2004`
    """
    # Comment before: this has to change. Now database should be a dictionary (cf Example) # TODO
    # Since 'keywords' is a dictionary, we cannot access it by index (e.g., database[0]).
    # Instead, we get the first value (file path) using list(keywords.values())[0].
    dir_path = list(keywords.values())[0]

    # Remove the filename portion to get the directory path.
    # This works by finding the last "/" in the file path and slicing up to that position.
    dir_path = dir_path[: dir_path.rindex("/") + 1]
    # Example: If dir_path = "path/to/data/v0.npy"
    # dir_path.rindex("/") returns 11 (the position of the last "/")
    # dir_path[: dir_path.rindex("/") + 1] gives "path/to/data/"

    try:
        if verbose >= 2:
            print("Loading iso...", end=" ")
        iso = np.load(dir_path + "iso.npy")
        if verbose >= 2:
            print("Done!")

        if verbose >= 2:
            print("Loading v0...", end=" ")
        v0 = np.load(dir_path + "v0.npy")
        if verbose >= 2:
            print("Done!")

        if verbose >= 2:
            print("Loading da...", end=" ")
        da = np.load(dir_path + "da.npy")
        if verbose >= 2:
            print("Done!")

        if verbose >= 2:
            print("Loading log_2gs...", end=" ")
        log_2gs = np.load(dir_path + "log_2gs.npy")
        if verbose >= 2:
            print("Done!")

        if verbose >= 2:
            print("Loading S0...", end=" ")
        S0 = np.load(dir_path + "S0.npy")
        if verbose >= 2:
            print("Done!")

        if verbose >= 2:
            print("Loading El...", end=" ")
        El = np.load(dir_path + "El.npy")
        if verbose >= 2:
            print("Done!")

        if verbose >= 2:
            print("Loading log_2vMm...", end=" ")
        log_2vMm = np.load(dir_path + "log_2vMm.npy")
        if verbose >= 2:
            print("Done!")

        if verbose >= 2:
            print("Loading na...", end=" ")
        na = np.load(dir_path + "na.npy")
        if verbose >= 2:
            print("Done!")

        # Create DataFrame
        df = pd.DataFrame(
            {
                "iso": iso,
                "wav": v0,
                "Pshft": da,
                "log_2gs": log_2gs,
                "Tdpair": na,
                "log_2vMm": log_2vMm,
                "int": S0,
                "El": El,
            }
        )

    except FileNotFoundError:
        raise FileNotFoundError("Could not find npy dataset in the given directory")

    return df
