"""
HDF5 to Parquet migration utility for RADIS spectroscopic databases.

Converts existing Vaex/PyTables HDF5 files in ~/.radisdb/ to sorted
Parquet files with zstd compression. Sorting by 'wav' (wavenumber)
enables predicate pushdown when using Polars.

Runs automatically on first use when engine='polars' is selected.
"""

import os
import logging
from pathlib import Path

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def convert_hdf5_to_parquet(
    hdf5_path,
    parquet_path=None,
    sort_by="wav",
    row_group_size=100_000,
    compression="zstd",
    validate=True,
    rtol=1e-14,
    verbose=True,
):
    """Convert a single HDF5 file to sorted Parquet.

    Parameters
    ----------
    hdf5_path : str or Path
        Path to source .hdf5 or .h5 file
    parquet_path : str or Path, optional
        Output path. If None, replaces .hdf5/.h5 extension with .parquet
    sort_by : str
        Column to sort by. Default 'wav' for wavenumber — this is important
        because Parquet stores min/max statistics per row group, and sorting
        by wavenumber means predicate pushdown can skip irrelevant chunks.
    row_group_size : int
        Number of rows per row group. 100k is a good balance between
        metadata overhead and filter granularity for HITRAN-sized data.
    compression : str
        Parquet compression. 'zstd' gives best ratio for spectroscopic data.
    validate : bool
        If True, read back and compare with original.
    rtol : float
        Relative tolerance for validation. 1e-14 is tight enough for
        spectroscopic accuracy (HITRAN stores 10 significant digits).
    verbose : bool
        Print progress messages.

    Returns
    -------
    str
        Path to the created Parquet file.
    """
    hdf5_path = Path(hdf5_path)
    if not hdf5_path.exists():
        raise FileNotFoundError(f"HDF5 file not found: {hdf5_path}")

    if parquet_path is None:
        parquet_path = hdf5_path.with_suffix(".parquet")
    parquet_path = Path(parquet_path)

    if parquet_path.exists():
        if verbose:
            logger.info(f"Parquet file already exists: {parquet_path}")
        return str(parquet_path)

    if verbose:
        logger.info(f"Converting {hdf5_path} -> {parquet_path}")

    # try vaex first (for vaex-format hdf5), fall back to pytables
    df = _read_hdf5(hdf5_path)

    if verbose:
        logger.info(f"  Loaded {len(df)} rows, {len(df.columns)} columns")

    # sort by wavenumber for optimal predicate pushdown
    if sort_by and sort_by in df.columns:
        df = df.sort_values(sort_by).reset_index(drop=True)
        if verbose:
            logger.info(f"  Sorted by '{sort_by}'")

    # write parquet
    df.to_parquet(
        parquet_path,
        engine="pyarrow",
        compression=compression,
        row_group_size=row_group_size,
        index=False,
    )

    file_size_mb = parquet_path.stat().st_size / 1e6
    if verbose:
        logger.info(f"  Written {file_size_mb:.1f} MB ({compression} compression)")

    # validate
    if validate:
        _validate_conversion(df, parquet_path, rtol, verbose)

    return str(parquet_path)


def _read_hdf5(path):
    """Read HDF5 file, trying vaex format first then pytables."""
    path = Path(path)

    # try vaex format first
    try:
        import vaex

        vdf = vaex.open(str(path))
        df = vdf.to_pandas_df()
        logger.debug(f"  Read as vaex HDF5: {path.name}")
        return df
    except Exception:
        pass

    # fall back to pytables (pandas)
    try:
        df = pd.read_hdf(path)
        logger.debug(f"  Read as pytables HDF5: {path.name}")
        return df
    except Exception:
        pass

    # last resort: h5py
    try:
        import h5py

        with h5py.File(path, "r") as f:
            data = {}
            # radis stores data in the 'df' group or root
            group = f["df"] if "df" in f else f
            if hasattr(group, "keys"):
                for key in group.keys():
                    if hasattr(group[key], "shape"):
                        data[key] = group[key][:]
            else:
                # try table format
                for key in f.keys():
                    if hasattr(f[key], "shape"):
                        data[key] = f[key][:]
        df = pd.DataFrame(data)
        logger.debug(f"  Read with h5py: {path.name}")
        return df
    except Exception as e:
        raise RuntimeError(
            f"Could not read HDF5 file {path}. "
            f"Tried vaex, pytables and h5py. Last error: {e}"
        )


def _validate_conversion(original_df, parquet_path, rtol, verbose):
    """Check that parquet file matches original data."""
    check_df = pd.read_parquet(parquet_path)

    if len(check_df) != len(original_df):
        raise RuntimeError(
            f"Row count mismatch: original {len(original_df)}, "
            f"parquet {len(check_df)}"
        )

    # check numeric columns
    for col in original_df.columns:
        if original_df[col].dtype in (np.float64, np.float32):
            orig = original_df[col].sort_values().reset_index(drop=True)
            conv = check_df[col].sort_values().reset_index(drop=True)
            if not np.allclose(orig, conv, rtol=rtol, equal_nan=True):
                raise RuntimeError(f"Data mismatch in column '{col}'")

    if verbose:
        logger.info("  Validation passed")


def migrate_radisdb(
    radisdb_path=None,
    sort_by="wav",
    verbose=True,
):
    """Convert all HDF5 files in ~/.radisdb/ to Parquet.

    Parameters
    ----------
    radisdb_path : str or Path, optional
        Path to radisdb directory. Default ~/.radisdb/
    sort_by : str
        Column to sort by. Default 'wav'.
    verbose : bool
        Print progress.

    Returns
    -------
    list of str
        Paths to created Parquet files.
    """
    if radisdb_path is None:
        radisdb_path = Path.home() / ".radisdb"

    radisdb_path = Path(radisdb_path)
    if not radisdb_path.exists():
        if verbose:
            logger.info(f"No radisdb directory found at {radisdb_path}")
        return []

    converted = []
    hdf5_files = list(radisdb_path.rglob("*.hdf5")) + list(
        radisdb_path.rglob("*.h5")
    )

    if verbose:
        logger.info(f"Found {len(hdf5_files)} HDF5 files in {radisdb_path}")

    for hdf5_file in hdf5_files:
        try:
            parquet_path = convert_hdf5_to_parquet(
                hdf5_file,
                sort_by=sort_by,
                verbose=verbose,
            )
            converted.append(parquet_path)
        except Exception as e:
            logger.warning(f"Failed to convert {hdf5_file}: {e}")

    if verbose:
        logger.info(f"Converted {len(converted)} / {len(hdf5_files)} files")

    return converted


def get_parquet_path(hdf5_path):
    """Get the corresponding .parquet path for a .hdf5 file.

    If parquet doesn't exist yet, converts it first.
    This is the main entry point — called from DataFileManager
    when engine='polars'.
    """
    hdf5_path = Path(hdf5_path)
    parquet_path = hdf5_path.with_suffix(".parquet")

    if not parquet_path.exists():
        convert_hdf5_to_parquet(hdf5_path, parquet_path)

    return str(parquet_path)
