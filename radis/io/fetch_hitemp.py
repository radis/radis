# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 22:40:51 2021

@author: erwan


https://stackoverflow.com/questions/55610891/numpy-load-from-io-bytesio-stream
https://stupidpythonideas.blogspot.com/2014/07/three-ways-to-read-files.html

"""

import os
import sys
from os.path import exists, expanduser, join
from time import time

import numpy as np
from numpy import DataSource

from radis.io.hitran import columns_2004
from radis.io.tools import _create_dtype, _get_linereturnformat, _ndarray2df

BASE_URL = "https://hitran.org/hitemp/data/bzip2format/"
SOURCE_FILES = {
    "CO2": "",  # NotImplemented
    "N2O": "04_HITEMP2019.par.bz2",
    "CO": "05_HITEMP2019.par.bz2",
    "CH4": "06_HITEMP2020.par.bz2",
    "NO": "08_HITEMP2019.par.bz2",
    "NO2": "10_HITEMP2019.par.bz2",
    "OH": "13_HITEMP2020.par.bz2",
}


def fetch_hitemp(molecule, local_databases="~/.radisdb/", verbose=True):
    """Stream HITEMP file from HITRAN website. Unzip and build a HDF5 file directly

    Parameters
    ----------
    molecule: `"CO2", "N2O", "CO", "CH4", "NO", "NO2", "OH"`
        HITEMP molecule
    local_databases: str
        where to download the files
    verbose: bool

    Returns
    -------
    path: str
        path of .h5 file where lines have been stored.

    """
    # TODO: add metadata (link, time-downloaded, min/max, etc.)
    # TODO: do not re-download if HDF5 exists  (currently DataSource is initialized)
    # TODO: clean DataSource downloads after HDF5 is generated

    local_databases = local_databases.replace("~", expanduser("~"))

    if molecule in ["H2O", "CO2"]:
        raise NotImplementedError(
            "Multiple files. Download manually on https://hitran.org/hitemp/ "
        )

    inputf = SOURCE_FILES[molecule]
    urlname = BASE_URL + inputf

    try:
        os.mkdir(local_databases)
    except OSError:
        pass
    else:
        if verbose:
            print("Created folder :", local_databases)

    output = join(local_databases, inputf.replace(".par.bz2", ".h5"))
    ds = DataSource(local_databases)

    if verbose:
        print(f"Downloading {inputf} for {molecule} in {output}")

    columns = columns_2004

    # Get linereturn (depends on OS, but file may also have been generated
    # on a different OS. Here we simply read the file to find out)
    with ds.open(urlname) as gfile:  # locally downloaded file

        dt = _create_dtype(columns, "a2")  # 'a2' allocates space to get \n or \n\r
        b = np.empty(161, dtype=dt)
        gfile.readinto(b)
        linereturnformat = _get_linereturnformat(b, columns)

    with ds.open(urlname) as gfile:  # locally downloaded file

        CHUNK = 1000  # TODO: adjust if needed (RAM dependant?)

        dt = _create_dtype(columns, linereturnformat)
        b = np.empty(
            dt.itemsize * CHUNK, dtype=dt
        )  # receives the HITRAN 160-format data.

        if exists(output):
            if verbose:
                print("Removing existing file ", output)
            os.remove(output)

        t0 = time()
        wmin = np.inf
        wmax = 0
        if verbose:
            print("Building", output)

        Nlines = 0
        for nbytes in iter(lambda: gfile.readinto(b), 0):

            if not b[-1]:
                # End of file within the chunk (but does not start with end of file)
                # so nbytes != 0
                b = get_last(b)

            df = _ndarray2df(b, columns, linereturnformat)

            df.to_hdf(
                output, "df", format="table", append=True, complib="blosc", complevel=9
            )

            wmin = np.min((wmin, df.wav.min()))
            wmax = np.max((wmax, df.wav.max()))
            Nlines += len(df)
            if verbose:
                sys.stdout.write(
                    "\r({0:.0f}s) Built {1} database from {2:.1f} to {3:.1f} cm-1 ({4} lines)".format(
                        time() - t0, molecule, wmin, wmax, Nlines
                    )
                )

            # Reinitialize for next read
            b = np.empty(
                dt.itemsize * CHUNK, dtype=dt
            )  # receives the HITRAN 160-format data.

    return output


#%%
def get_last(b):
    """ Get non-empty lines of a chunk b, parsing the bytes """

    element_length = np.vectorize(lambda x: len(x.__str__()))(b)
    non_zero = element_length > element_length[-1]
    threshold = non_zero.argmin() - 1
    assert (non_zero[: threshold + 1] == 1).all()
    assert (non_zero[threshold + 1 :] == 0).all()
    return b[non_zero]


#%%

if __name__ == "__main__":

    fetch_hitemp("OH")
