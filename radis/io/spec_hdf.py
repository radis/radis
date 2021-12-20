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
    """Convert Spectrum s to HDF5 file"""

    df_arrays = s.to_pandas(copy=False)

    # Store
    mgr = HDF5Manager(engine=engine)

    # Store spectral arrays and lines as Datasets

    mgr.write(file, df_arrays, key="arrays", append=False, data_columns=["wavespace"])
    mgr.add_metadata(file, s.units, key="arrays")

    if s.lines:
        mgr.write(file, s.lines, key="lines", append=True)

    # Store the rest as metadata

    if s.populations:
        mgr.add_metadata(
            file, s.populations, key="populations", create_empty_dataset=True
        )

    mgr.add_metadata(file, s.conditions, key="conditions", create_empty_dataset=True)

    if s.references:
        mgr.add_metadata(
            file, s.references, key="references", create_empty_dataset=True
        )


def hdf2spec(file, wmin=None, wmax=None, wunit=None, columns=None, engine="pytables"):
    """Read HDF5 file into a Spectrum object"""
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
            #  Selection is done after opening the file in vaex
            where = None
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
        if wmin is not None and wmax is not None:
            df2.select((df2.wavespace > wmin) & (df2.wavespace < wmax))
            selection = True
        elif wmin is not None:
            df2.select(df2.wavespace > wmin)
            selection = True
        elif wmax is not None:
            df2.select(df2.wavespace < wmax)
            selection = True
        else:
            selection = False
        df2 = df2.to_pandas_df(selection=selection, column_names=columns)

    units2 = mgr.read_metadata(file, key="arrays")
    conditions2 = mgr.read_metadata(file, key="conditions")
    try:
        populations2 = mgr.read_metadata(file, key="populations")
    except KeyError:
        populations2 = None
    try:
        references2 = mgr.read_metadata(file, key="references")
    except KeyError:
        references2 = None

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
