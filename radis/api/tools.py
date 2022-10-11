# -*- coding: utf-8 -*-
"""Created on Fri Jul  6 13:52:04 2018.

@author: erwan
"""
# TODO refactor : rename this file as hitran_utils.py

import numpy as np
import pandas as pd


def parse_hitran_file(fname, columns, count=-1):
    """Parse a file under HITRAN ``par`` format. Parsing is done in binary
    format with :py:func:`numpy.fromfile` so it's as fast as possible.

    Parameters
    ----------
    fname: str
        filename.
    columns: dict
        list of columns and their format.

    Other Parameters
    ----------------
    count: int
        number of lines to read. If ``-1`` reads all file.

    Returns
    -------
    df: pandas DataFrame
        dataframe with lines.

    See Also
    --------
    Used in :py:func:`~radis.api.hitranapi.hit2df` and :py:func:`~radis.api.cdsdapi.cdsd2df`
    """

    # To be faster, we read file totally in bytes mode with fromfiles. But that
    # requires to properly decode the line return character:

    # We cannot simply infer it from the OS :
    # problem arise when file was written in an OS and read in another OS (for instance,
    # line return characters are not converted when read from .egg files). Here
    # we read the first line and infer the line return character for it
    data = _read_hitran_file(
        fname, columns, count=1, linereturnformat="a2"
    )  # 'a2' allocates space to get \n or \n\r linereturn formats
    linereturnformat = _get_linereturnformat(data, columns, fname)

    # Now re-read with correct line return character
    data = _read_hitran_file(fname, columns, count, linereturnformat)

    # Return a Pandas dataframe
    return _ndarray2df(data, columns, linereturnformat)


def _get_linereturnformat(data, columns, fname=""):
    """
    Get line return character & format (size).

    Notes
    -----

    We cannot simply infer it from the OS :
    problem arise when file was written in an OS and read in another OS (for instance,
    line return characters are not converted when read from .egg files). Here
    we read the first line and infer the line return character for it
    """
    # fname just for the error message

    # get format (size) of line return
    from radis.misc.basics import to_str

    linereturn = to_str(data[0][-1])
    if to_str("\r\n") in linereturn:
        linereturnformat = "a2"
    elif to_str("\n") in linereturn or to_str("\r") in linereturn:
        linereturnformat = "a1"
    else:
        raise ValueError(
            "Unknown Line return format: {0}. Check that your file {1} has the HITRAN format. First line : {2}".format(
                linereturn, fname, data[0]
            )
        )

    return linereturnformat


def _ndarray2df(data, columns, linereturnformat):
    """ """

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


def _read_hitran_file(fname, columns, count, linereturnformat):
    """
    Returns
    -------

    data: np.ndarray
    """

    dt = _create_dtype(columns, linereturnformat)

    return np.fromfile(fname, dtype=dt, count=count)


def _create_dtype(columns, linereturnformat):
    """Create dtype from specific columns.

    Returns
    -------
    dt: dtypes
    """
    # ... Create a dtype with the binary data format and the desired column names
    dtype = [(k, c[0]) for (k, c) in columns.items()] + [
        ("_linereturn", linereturnformat)
    ]
    # ... _linereturn is to capture the line return symbol. We delete it afterwards
    return _format_dtype(dtype)


def _format_dtype(dtype):
    """Format dtype from specific columns.

    Crash with hopefully helping error message
    """
    try:
        dt = np.dtype([(str(k), c) for k, c in dtype])
        # Note: dtype names cannot be `unicode` in Python2. Hence the str()
    except TypeError as err:
        # Cant read database. Try to be more explicit for user before crashing
        try:
            print("Data type")
            print("-" * 30)
            for (k, c) in dtype:
                print(str(k), "\t\t", c)
            print("-" * 30)
        finally:
            raise TypeError("Couldnt read datatype. See details above.") from err
    return dt


def _cast_to_dtype(data, dtype):
    """Cast array to certain type, crash with hopefull helping error message.

    Return casted data.

    Parameters
    ----------
    data: array to cast
    dtype: (ordered) list of (param, type)
    """

    dt = _format_dtype(dtype)
    try:
        data = np.array(data, dtype=dt)
    except ValueError as err:
        try:
            # Cant read database. Try to be more explicit for user before crashing
            # ... identify faulty row
            print("Cant cast data to specific dtype. Trying row by row:")
            for r in range(len(data)):
                try:
                    np.array(data[r], dtype=dt)
                except ValueError:
                    break
            print(f"Error may be on row {r}:")
            print("-" * 30)
            for i in range(len(data[r])):
                print(i, dtype[i], "\t\t", np.array(data[r][i], dtype=dt[i]))
            print("-" * 30)
            print(">>> Next param:", dtype[i], ". Value:", data[r][i], "\n")
        finally:
            raise ValueError(
                "Cant cast data to specific dtype. Tried column by column. See results above"
            ) from err

    return data


def drop_object_format_columns(df, verbose=True):
    """Remove 'object' columns in a pandas DataFrame.

    They are not useful to us at this time, and they slow down all
    operations (as they are converted to 'object' in pandas DataFrame).
    If you want to keep them, better convert them to some numeric values
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
    """Return P, Q, R in column ``branch`` with -1, 0, 1 to get a fully numeric
    database. This improves performances quite a lot, as Pandas doesnt have a
    fixed-string dtype hence would use the slow ``object`` dtype.

    Parameters
    ----------

    df: pandas Dataframe
        ``branch`` must be a column name.


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
