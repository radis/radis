#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Wrapper to fetch line database from HITRAN with Astroquery [1]_
(based on [HAPI]_)

.. note::
    if using, cite [HAPI]_ and [HITRAN-2020]_


References
----------

.. [R1] `Astroquery <https://astroquery.readthedocs.io>`_


-------------------------------------------------------------------------------

"""

import os
import sys
from os.path import exists, isfile, join

import numpy as np
import pandas as pd

from radis.db.classes import get_molecule, get_molecule_identifier

try:
    from .cache_files import check_cache_file, get_cache_file, save_to_hdf
except ImportError:
    from radis.api.cache_files import check_cache_file, get_cache_file, save_to_hdf
from radis.misc import is_float
from radis.misc.printer import printr

CACHE_FILE_NAME = "tempfile_{molecule}_{isotope}_{wmin:.2f}_{wmax:.2f}.h5"


def fetch_astroquery(
    molecule,
    isotope,
    wmin,
    wmax,
    verbose=True,
    cache=True,
    expected_metadata={},
    engine="pytables-fixed",
):
    """Download a HITRAN line database to a Pandas DataFrame.

    Wrapper to the fetch function of Astroquery [1]_ (itself based on [HAPI]_)

    .. note::
        if using, cite [HAPI]_ and [HITRAN-2020]_

    Parameters
    ----------
    molecule: str, or int
        molecule name or identifier
    isotope: int
        isotope number
    wmin, wmax: float  (cm-1)
        wavenumber min and max

    Other Parameters
    ----------------
    verbose: boolean
        Default ``True``
    cache: boolean or ``'regen'``
        if ``True``, tries to find a ``.h5`` cache file in the Astroquery
        :py:attr:`~astroquery.query.BaseQuery.cache_location`, that would match
        the requirements. If not found, downloads it and saves the line dataframe
        as a ``.h5`` file in the Astroquery.
        If ``'regen'``, delete existing cache file to regerenate it.
    expected_metadata: dict
        if ``cache=True``, check that the metadata in the cache file correspond
        to these attributes. Arguments ``molecule``, ``isotope``, ``wmin``, ``wmax``
        are already added by default.

    References
    ----------
    .. [1] `Astroquery <https://astroquery.readthedocs.io>`_

    See Also
    --------
    :py:meth:`astroquery.hitran.core.Hitran.query_lines_async`,
    :py:attr:`astroquery.query.BaseQuery.cache_location`

    """
    from astropy import units as u
    from astroquery.hitran import Hitran

    # Check input
    if not is_float(molecule):
        mol_id = get_molecule_identifier(molecule)
    else:
        mol_id = molecule
        molecule = get_molecule(mol_id)
    assert is_float(isotope)

    empty_range = False

    if cache:
        # Cache file location in Astroquery cache
        # TODO: move full HITRAN databases in ~/radisdb cache like io/hitemp/fetch_hitemp ?
        fcache = join(
            Hitran.cache_location,
            CACHE_FILE_NAME.format(
                **{"molecule": molecule, "isotope": isotope, "wmin": wmin, "wmax": wmax}
            ),
        )
        # ... Update metadata with physical properties from the database.
        expected_metadata.update(
            {"molecule": molecule, "isotope": isotope, "wmin": wmin, "wmax": wmax}
        )
        if cache == "regen":
            if exists(fcache):
                if verbose:
                    print(f"Cache file {fcache} deleted to be regenerated")
                os.remove(fcache)
        else:
            # Load cache file if valid
            check_cache_file(
                fcache=fcache,
                use_cached=cache,
                expected_metadata=expected_metadata,
                compare_as_close=["wmin", "wmax"],
                verbose=verbose,
                engine=engine,
            )
            if exists(fcache):
                try:
                    return get_cache_file(fcache, verbose=verbose, engine=engine)
                except Exception as err:
                    if verbose:
                        printr(
                            "Problem reading cache file {0}:\n{1}\nDeleting it!".format(
                                fcache, str(err)
                            )
                        )
                    os.remove(fcache)

    # Download using the astroquery library
    try:
        response = Hitran.query_lines_async(
            molecule_number=mol_id,
            isotopologue_number=isotope,
            min_frequency=wmin / u.cm,
            max_frequency=wmax / u.cm,
        )
    except KeyError as err:
        raise KeyError(
            str(err)
            + " <<w this error occured in Astroquery. Maybe these molecule "
            + "({0}) and isotope ({1}) are not supported".format(molecule, isotope)
        ) from err

    # Deal with usual errors
    if response.status_code == 404:
        # Maybe there are just no lines for this species in this range
        # In that case we usually end up with errors like:

        # (<class 'Exception'>, Exception('Query failed: 404 Client Error:
        # Not Found for url: http://hitran.org/lbl/api?numax=25000&numin=19000&iso_ids_list=69\n',),
        # <traceback object at 0x7f0967c91708>)

        if response.reason == "Not Found":
            # Let's bet it's just that there are no lines in this range
            empty_range = True
            if verbose:
                print(
                    (
                        "No lines for {0} (id={1}), iso={2} in range {3:.2f}-{4:.2f}cm-1. ".format(
                            molecule, mol_id, isotope, wmin, wmax
                        )
                    )
                )
        else:
            raise ValueError(
                "An error occured during the download of HITRAN files "
                + "for {0} (id={1}), iso={2} between {3:.2f}-{4:.2f}cm-1. ".format(
                    molecule, mol_id, isotope, wmin, wmax
                )
                + "Are you online?\n"
                + "See details of the error below:\n\n {0}".format(response.reason)
            )
    elif response.status_code == 500:

        raise ValueError(
            "{0} while querying the HITRAN server: ".format(response.status_code)
            + "\n\n{0}".format(response.text)
        )

    # Process response

    # Rename columns from Astroquery to RADIS format
    # TODO : define RADIS database format somewhere else; with description of the column names.
    rename_columns = {
        "molec_id": "id",
        "local_iso_id": "iso",
        "nu": "wav",
        "sw": "int",
        "a": "A",
        "gamma_air": "airbrd",
        "gamma_self": "selbrd",
        "elower": "El",
        "n_air": "Tdpair",
        "delta_air": "Pshft",
        "global_upper_quanta": "globu",
        "global_lower_quanta": "globl",
        "local_upper_quanta": "locu",
        "local_lower_quanta": "locl",
        "line_mixing_flag": "lmix",
        "gp": "gp",
        "gpp": "gpp",
    }

    if not empty_range:
        try:
            tbl = Hitran._parse_result(response)
        except ValueError as err:
            raise ValueError(
                "Error while parsing HITRAN output : {}".format(response.text)
            ) from err
        df = tbl.to_pandas()
        df = df.rename(columns=rename_columns)
    else:
        df = pd.DataFrame(columns=list(rename_columns.values()))

    # Cast type to float64
    cast_type = {
        "wav": np.float64,
        "int": np.float64,
        "A": np.float64,
        "airbrd": np.float64,
        "selbrd": np.float64,
        "El": np.float64,
        "Tdpair": np.float64,
        "Pshft": np.float64,
    }
    for c, typ in cast_type.items():
        df[c] = df[c].astype(typ)

    # cached file mode but cached file doesn't exist yet (else we had returned)
    if cache:
        new_metadata = {
            "molecule": molecule,
            "isotope": isotope,
            "wmin": wmin,
            "wmax": wmax,
        }
        if verbose:
            print(
                "Generating cache file {0} with metadata :\n{1}".format(
                    fcache, new_metadata
                )
            )
        from radis import __version__

        try:
            save_to_hdf(
                df,
                fcache,
                metadata=new_metadata,
                version=__version__,
                key="df",
                overwrite=True,
                verbose=verbose,
                engine=engine,
            )
        except PermissionError:
            if verbose:
                print(sys.exc_info())
                print("An error occured in cache file generation. Lookup access rights")
            pass

    return df


def _fix_astroquery_file_format(filename):
    """
    Notes
    -----

    On some OS the astroquery lookup function may add extra lines. See:
    https://github.com/astropy/astroquery/issues/1189

    In the meantime, we discard all empty lines here.

    """

    if not isfile(filename):
        raise "{} does not exist "
    with open(filename) as filehandle:
        lines = filehandle.readlines()
        non_empty_lines = [
            l for l in lines if len(l) > 2
        ]  # > 2 because there may be some line return characters

    if len(lines) != len(non_empty_lines):
        # Re-write file
        with open(filename, "w") as filehandle:
            filehandle.writelines(non_empty_lines)


if __name__ == "__main__":
    from radis.test.io.test_query import _run_testcases

    _run_testcases(verbose=True)
