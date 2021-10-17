# -*- coding: utf-8 -*-
"""
Functions to read/write Spectrum objects under HDF5 format
"""

try:
    from .hdf5 import HDF5Manager
except ImportError:
    from radis.io.hdf5 import HDF5Manager


# Convert


def spec2hdf(s, file, engine="pytables"):
    """ Convert Spectrum s to HDF5 file"""

    df_arrays = s.to_pandas(copy=False)

    # Store
    mgr = HDF5Manager(engine=engine)

    # Store spectral arrays and lines as Datasets

    mgr.write(file, df_arrays, key="arrays", append=False, data_columns=["wavespace"])

    if s.lines:
        mgr.write(file, s.lines, key="lines", append=True)

    # Store the rest as metadata

    metadata = {}
    metadata["units"] = s.units

    if s.populations:
        metadata["populations"] = s.populations

    metadata["conditions"] = s.conditions

    if s.references:
        metadata["references"] = s.references

    mgr.add_metadata(file, metadata, key="arrays")


def hdf2spec(file, wmin=None, wmax=None, wunit=None, columns=None, engine="pytables"):
    """ Read HDF5 file into a Spectrum object"""
    from radis import Spectrum

    mgr = HDF5Manager(engine=engine)

    # Range selection
    # TODO : convert if stored in wavelength/wavenumber
    if wmin is not None or wmax is not None:
        if engine == "pytables":
            where = []
            if wmin is not None:
                where.append(f"wavespace > {wmin}")
            if wmax is not None:
                where.append(f"wavespace < {wmax}")
        elif engine == "vaex":
            #  Selection is done after opening the file time in vaex
            pass
        else:
            raise NotImplementedError(engine)
    else:
        where = None

    # Column selection
    # ... always load 'wavespace'
    if columns and "wavespace" not in columns:
        columns.append("wavespace")

    df2 = mgr.load(file, columns=columns, where=where, key="arrays")

    # Convert to Pandas df (if vaex)
    if engine == "vaex":
        if wmin is not None or wmax is not None:
            raise NotImplementedError(engine)
        df2 = df2.to_pandas_df()

    metadata = mgr.read_metadata(file, key="arrays")

    units2 = metadata["units"]
    conditions2 = metadata["conditions"]
    populations2 = metadata.get("populations", None)
    references2 = metadata.get("references", None)

    try:
        lines2 = mgr.load(file, key="lines")
    except (KeyError, OSError):  # KeyError in PyTables ; OSError in Vaex
        lines2 = None

    s2 = Spectrum(
        {k: v.values for k, v in dict(df2).items()},
        lines=lines2,
        units=units2,
        conditions=conditions2,
        references=references2,
        populations=populations2,
    )

    return s2
