# -*- coding: utf-8 -*-
"""
Created on Fri Jul  6 13:52:04 2018

@author: erwan
"""

from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import pandas as pd
from six.moves import range
from six.moves import zip


def parse_hitran_file(fname, columns, count):
    """ Parse a file under HITRAN ``par`` format. Parsing is done in binary 
    format with :py:func:`numpy.fromfile` so it's as fast as possible.

    Parameters
    ----------

    fname: str
        filename

    columns: dict
        list of columns and their format

    count: int
        number of lines to read

    Returns
    -------

    df: pandas DataFrame
        dataframe with lines

    Notes
    -----

    Part common to hit2df and cdsd2df

    """

    # To be faster, we read file totally in bytes mode with fromfiles. But that
    # requires to properly decode the line return character:

    # problem arise when file was written in an OS and read in another OS (for instance,
    # line return characters are not converted when read from .egg files). Here
    # we read the first line and infer the line return character for it

    # ... Create a dtype with the binary data format and the desired column names
    dtype = [(k, c[0]) for (k, c) in columns.items()] + [("_linereturn", "a2")]
    # ... _linereturn is to capture the line return symbol. We delete it afterwards
    dt = _format_dtype(dtype)
    data = np.fromfile(fname, dtype=dt, count=1)  # just read the first line

    # get format of line return
    from radis.misc.basics import to_str

    linereturn = to_str(data[0][-1])
    if to_str("\r\n") in linereturn:
        linereturnformat = "a2"
    elif to_str("\n") in linereturn or to_str("\r") in linereturn:
        linereturnformat = "a1"
    else:
        raise ValueError(
            "Unknown `Line return` format: {0}. Check that your file {0} has the HITRAN format.".format(
                linereturn, fname
            )
        )

    # Now re-read with correct line return character

    # ... Create a dtype with the binary data format and the desired column names
    dtype = [(k, c[0]) for (k, c) in columns.items()] + [
        ("_linereturn", linereturnformat)
    ]
    # ... _linereturn is to capture the line return symbol. We delete it afterwards
    dt = _format_dtype(dtype)
    data = np.fromfile(fname, dtype=dt, count=count)

    # ... Cast to new type
    # This requires to recast all the data already read, but is still the fastest
    # method I found to read a file directly (for performance benchmark see
    # CDSD-HITEMP parser)
    newtype = [c[0] if (c[1] == str) else c[1] for c in columns.values()]
    dtype = list(zip(list(columns.keys()), newtype)) + [
        ("_linereturn", linereturnformat)
    ]
    data = _cast_to_dtype(data, dtype)

    # %% Create dataframe
    df = pd.DataFrame(data.tolist(), columns=list(columns.keys()) + ["_linereturn"])

    # Delete dummy column than handled the line return character
    del df["_linereturn"]

    # Update format
    for k, c in columns.items():
        if c[1] == str:
            df[k] = df[k].str.decode("utf-8")

    # Strip whitespaces around PQR columns (due to 2 columns jumped)
    if "branch" in df:  # (only in CDSD)
        df["branch"] = df.branch.str.strip()

    return df


def _format_dtype(dtype):
    """ Format dtype from specific columns. Crash with hopefully helping error message """

    try:
        dt = np.dtype([(str(k), c) for k, c in dtype])
        # Note: dtype names cannot be `unicode` in Python2. Hence the str()
    except TypeError:
        # Cant read database. Try to be more explicit for user
        print("Data type")
        print(("-" * 30))
        for (k, c) in dtype:
            print((str(k), "\t", c))
        print(("-" * 30))
        raise
    return dt


def _cast_to_dtype(data, dtype):
    """ Cast array to certain type, crash with hopefull helping error message.
    Return casted data


    Parameters    
    ----------

    data: array to cast

    dtype: (ordered) list of (param, type)

    """

    dt = _format_dtype(dtype)

    try:
        data = np.array(data, dtype=dt)
    except ValueError:
        try:
            # Cant read database. Try to be more explicit for user
            print("Cant cast data to specific dtype. Trying column by column:")
            print(("-" * 30))
            for i in range(len(data[0])):
                print((dtype[i], "\t", np.array(data[0][i], dtype=dt[i])))
            print(("-" * 30))
        except ValueError:
            print((">>> Next param:", dtype[i], ". Value:", data[0][i], "\n"))
            raise ValueError(
                "Cant cast data to specific dtype. Tried column by column. See results above"
            )

    return data


def drop_object_format_columns(df, verbose=True):
    """ Remove 'object' columns in a pandas DataFrame. They are not useful to us at this 
    time, and they slow down all operations (as they are converted to 'object'
    in pandas DataFrame). If you want to keep them, better convert them to 
    some numeric values
    """

    objects = [k for k, v in df.dtypes.items() if v == object]
    for k in objects:
        del df[k]
    if verbose >= 2 and len(objects) > 0:
        print(
            (
                "The following columns had the `object` format and were removed: {0}".format(
                    objects
                )
            )
        )
    return df


def replace_PQR_with_m101(df):
    """ Return P, Q, R in column ``branch`` with -1, 0, 1 to get a fully numeric 
    database. This improves performances quite a lot, as Pandas doesnt have a 
    fixed-string dtype hence would use the slow ``object`` dtype. 
    
    Parameters
    ----------
    
    df: pandas Dataframe
        ``branch`` must be a column name
    
    Returns
    -------
    
    None:
        ``df`` is is modified in place 
    """

    # Note: somehow pandas updates dtype automatically. We have to check
    # We just have to replace the column:
    if df.dtypes["branch"] != np.int64:
        new_col = df["branch"].replace(["P", "Q", "R"], [-1, 0, 1])
        df["branch"] = new_col

    assert df.dtypes["branch"] == np.int64
