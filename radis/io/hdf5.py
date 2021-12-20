# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 21:27:15 2021

@author: erwan
"""

import os
import sys
from os.path import exists, expanduser, splitext
from time import time

import h5py
import pandas as pd
from tables.exceptions import NoSuchNodeError


def update_pytables_to_vaex(fname, remove_initial=False, verbose=True, key="df"):
    """Convert a HDF5 file generated from PyTables to a
    Vaex-friendly HDF5 format, preserving metadata"""
    import vaex

    if fname.endswith(".h5"):
        fname_vaex = fname.replace(".h5", ".hdf5")
    else:
        fname_vaex = fname

    if verbose:
        print(f"Auto-updating {fname} to a Vaex-compatible HDF5 file {fname_vaex}")
    df = pd.read_hdf(fname)
    df = vaex.from_pandas(df)

    pytables_manager = HDF5Manager(engine="pytables")

    # Read metadata
    file_metadata = pytables_manager.read_metadata(fname)

    # Write Vaex file
    df.export_hdf5(fname_vaex)
    df.close()  # try to fix file not closed()  TODO: remove?
    del df  # same TODO

    vaex_manager = HDF5Manager(engine="vaex")
    vaex_manager.add_metadata(fname_vaex, file_metadata)

    if verbose:
        print(f"Converted to Vaex's HDF5 format {fname_vaex}")

    if remove_initial and fname != fname_vaex:
        # Remove initial file
        os.remove(fname)
        if verbose:
            print(f"Deleting {fname}")

    return fname_vaex


class HDF5Manager(object):
    def __init__(self, engine=None):
        """Class to handle all memory-mapping-librairies with one common API

        All functions may not be fully implemetned, will raise a NotImplementedError
        if that's not the case.

        Librairies ::

            'vaex'     > HDF5,  column-based
            'pytables' > Pandas's HDF5,  row-based
            'h5py'     > HDF5

        Functions ::

            add_metadata
            read_metadata
            write
            load
            guess_engine

        Examples
        --------
        ::

            file = 'CO.hdf5'
            from radis.io.hdf5 import HDF5Manager
            engine = HDF5Manager.guess_engine(file)
            mgr = HDF5Manager(engine)
            mgr.read_metadata(file)


        """
        self.engine = engine
        self._temp_batch_files = (
            []
        )  # list of batch files when writing by part in vaex mode

    def open(self, file, mode="w"):
        if self.engine == "pytables":
            return pd.HDFStore(file, mode=mode, complib="blosc", complevel=9)
        elif self.engine == "vaex":
            import vaex

            return vaex.open(file)
        else:
            raise NotImplementedError(self.engine)

    def write(
        self,
        file,
        df,
        append=True,
        key="default",
        format="table",
        data_columns=["iso", "wav"],
    ):
        """Write dataframe ``df`` to ``file``

        Other Parameters
        ----------------
        key: str
            group to write to. If ``None``, write at root level. If ``'default'``,
            use engine's default (`/table` for `'vaex'`, `df` for `pytables`,
            root for `h5py` )
        data_columns : list
            only these column names will be searchable directly on disk to
            load certain lines only. See :py:func:`~radis.io.hdf5.hdf2df`
        """
        file = expanduser(file)
        if self.engine == "pytables":
            if key == "default":
                key = "df"
            with self.open(file, "a" if append else "w") as f:
                f.put(
                    key=key,
                    value=df,
                    append=append,
                    format=format,
                    data_columns=data_columns,
                )
        elif self.engine == "pytables-fixed":
            assert not append
            # export dataframe
            df.to_hdf(file, key, format="fixed", mode="w", complevel=9, complib="blosc")
        elif self.engine == "vaex":
            if key == "default":
                key = r"/table"
            import vaex

            if append == True:
                # In vaex we cannot append. Here we write lots of small files then combine them.
                # self.combine_temp_batch_files() should be combined at the end.
                base, ext = splitext(file)
                i = 0
                temp_batch_file = base + "_temp" + str(i).zfill(6) + ext
                while temp_batch_file in self._temp_batch_files:
                    i += 1
                    temp_batch_file = base + "_temp" + str(i).zfill(6) + ext
                file = temp_batch_file
                # Check no remaining one from a non-cleaned previous run:
                if exists(file):
                    from radis.misc.printer import printr

                    printr(f"Temp file {temp_batch_file} already exists: deleting it")
                    os.remove(file)
                self._temp_batch_files.append(file)
            # Write:
            try:
                df.export_hdf5(file, group=key, mode="w")
            except AttributeError:  # case where df is not a Vaex dataFrame but (likely) a Pandas Dataframe
                vaex.from_pandas(df).export_hdf5(file, group=key, mode="w")
        else:
            raise NotImplementedError(self.engine)
            # h5py is not designed to write Pandas DataFrames

    def combine_temp_batch_files(self, file, key="default", sort_values=None):
        """Combine all batch files in ``self._temp_batch_files`` into one.
        Removes all batch files.
        """
        file = expanduser(file)
        if self.engine == "vaex":
            if len(self._temp_batch_files) == 0:
                raise ValueError(f"No batch temp files were written for {file}")
            if key == "default":
                key = r"/table"
            import vaex

            df = vaex.open(self._temp_batch_files, group=key)
            if sort_values:
                df.sort(by=sort_values).export_hdf5(file, group=key, mode="w")
            else:
                df.export_hdf5(file, group=key, mode="w")
            df.close()
        self._close_temp_batch_files()

    def _close_temp_batch_files(self):
        for i in range(len(self._temp_batch_files) - 1, -1, -1):
            os.remove(self._temp_batch_files[i])
            del self._temp_batch_files[i]

    def __del__(self):
        """ clean before deleting"""
        if len(self._temp_batch_files) > 0:
            from radis.misc.printer import printr

            printr(
                "Warning : {0} files were written but not combined. Deleting them.\n\n{1}".format(
                    len(self._temp_batch_files), self._temp_batch_files
                )
            )
            self._close_temp_batch_files()

    def load(
        self,
        fname,
        columns=None,
        where=None,
        key="default",
        none_if_empty=False,
        **store_kwargs,
    ):
        """
        Parameters
        ----------
        columns: list of str
            list of columns to load. If ``None``, returns all columns in the file.
        where: list of str
            filtering conditions. Ex::

                "wav > 2300"

        Other Parameters
        ----------------
        key: str
            group to load from. If ``None``, load from root level. If ``'default'``,
            use engine's default (`/table` for `'vaex'`, `df` for `pytables`,
            root for `h5py` )

        Returns
        -------
        pd.DataFrame or vaex.DataFrame
        """

        if self.engine in ["pytables", "pytables-fixed"]:
            fname = expanduser(fname)
            if key == "default":
                key = "df"
            try:
                df = pd.read_hdf(
                    fname, columns=columns, where=where, key=key, **store_kwargs
                )
            except TypeError as err:
                if "reading from a Fixed format store" in str(err):
                    raise TypeError(
                        f"radis.io.hdf5.hdf2df can only be used to load specific HDF5 files generated in a 'Table' which allows to select only certain columns or rows. Here the file {fname} is in 'Fixed' format. Regenerate it ? If it's a cache file of a .par file, load the .par file directly ?"
                    )
                elif "cannot create a storer if the object is not existing" in str(err):
                    raise TypeError(
                        f"Missing group `{key}` in {fname}. Maybe the file has been generated by a different HDF5 library than Pytables. Try using `engine='vaex'` in the calling function (hdf2df, etc.)"
                    )
                else:
                    raise
            except NoSuchNodeError as err:
                # Probably tried to read a Vaex/h5py HDF5 file forcing "engine='pytables'"
                raise AttributeError(
                    f"file {fname} does not seem to have been generated by Pytables. Try using `engine='vaex'` in the calling function (hdf2df, etc.)"
                ) from err

        elif self.engine == "vaex":
            if key == "default":
                key = r"/table"

            import vaex

            # Open file
            assert len(store_kwargs) == 0
            fname_list = fname if isinstance(fname, list) else [fname]
            fname_list = [expanduser(f) for f in fname_list]
            # vaex can open several files at the same time
            # First check group exists
            for fname in fname_list:
                try:
                    with h5py.File(fname, "r") as f:
                        if key and key not in f:
                            raise KeyError(
                                key
                            )  # before vaex raises the same error; which prints things on the console that cannot be caught
                except (FileNotFoundError, OSError) as err:
                    # error message with suggestion on how to convert from existing file
                    for f in fname_list:
                        if exists(f.replace(".hdf5", ".h5")):
                            raise FileNotFoundError(
                                f"`{f}` not found but `{f.replace('.hdf5', '.h5')}` exists (probably a row-based pytables HDF5 file). Try (1) using engine='pytables' in the calling function (`hdf2df`, `fetch_hitemp`, etc.)  ; (2) delete the file to re-download and re-parse it (this may take a lot of time !) ;  or (3, recommended) set `import radis; radis.config['AUTO_UPDATE_DATABASE']= True` in your script to auto-update to Vaex HDF5 file"
                            ) from err
                    raise

            # Now, open with vaex
            try:
                df = vaex.open(fname_list, group=key)
            except OSError as err:
                raise OSError(
                    f"Cannot read {fname}, group `{key}` with Vaex HDF5 library (column-based). It may be a file generated by pytables (row-based). Try (1) using engine='pytables' in the calling function (`hdf2df`, `fetch_hitemp`, etc.)  ; (2) delete the file to re-download and re-parse it (this may take a lot of time !) ;  or (3, recommended) set `import radis; radis.config['AUTO_UPDATE_DATABASE'] = True` in your script to auto-update to Vaex HDF5 file"
                ) from err

            return df

        elif self.engine == "h5py":
            fname = expanduser(fname)
            # TODO: define default key ?
            if key == "default":
                key = None

            with h5py.File(fname, "r") as f:
                if key is None:  # load from root level
                    load_from = f
                else:
                    load_from = f[key]
                out = {}
                for k in load_from.keys():
                    out[k] = f[k][()]
            return pd.DataFrame(out)

        else:
            raise NotImplementedError(self.engine)

        return df

    def add_metadata(
        self, fname: str, metadata: dict, key="default", create_empty_dataset=False
    ):
        """
        Parameters
        ----------
        fname: str
            filename
        metadata: dict
            dictionary of metadata to add in group ``key``
        key: str
            group to add metadata to. If ``None``, add at root level. If ``'default'``,
            use engine's default (`/table` for `'vaex'`, `df` for `pytables`,
            root for `h5py` )

        Other Parameters
        ----------------
        create_empty_dataset: bool
            if True, create an empty dataset to store the metadata as attribute

        """
        from radis.io.cache_files import _h5_compatible

        if self.engine in ["pytables", "pytables-fixed"]:
            fname = expanduser(fname)
            if key == "default":
                key = "df"
            with pd.HDFStore(fname, mode="a", complib="blosc", complevel=9) as f:
                if create_empty_dataset:
                    assert self.engine != "pytables-fixed"
                    # Not possible to create an empty node directly, so we create a dummy group  # TODO ?
                    f.put(key, pd.Series([]))
                f.get_storer(key).attrs.metadata = metadata

        elif self.engine == "h5py":
            fname = expanduser(fname)
            if key == "default":
                key = None
            with h5py.File(fname, "a") as hf:
                if create_empty_dataset:
                    assert key is not None
                    hf.create_dataset(key, dtype="f")
                if key is None:  # add metadta at root level
                    hf.attrs.update(_h5_compatible(metadata))
                else:
                    hf[key].attrs.update(_h5_compatible(metadata))

        elif self.engine == "vaex":
            if key == "default":
                key = r"/table"
            # Should be able to deal with multiple files at a time
            if isinstance(fname, list):
                assert isinstance(metadata, list)
                for f, m in zip(fname, metadata):
                    f = expanduser(f)
                    with h5py.File(f, "a") as hf:
                        if create_empty_dataset:
                            assert key is not None
                            hf.create_dataset(key, dtype="f")
                        if key is None:  # add metadta at root level
                            hf.attrs.update(_h5_compatible(m))
                        else:
                            hf[key].attrs.update(_h5_compatible(m))
            else:
                fname = expanduser(fname)
                with h5py.File(fname, "a") as hf:
                    if create_empty_dataset:
                        assert key is not None
                        hf.create_dataset(key, dtype="f")
                    if key is None:  # add metadta at root level
                        hf.attrs.update(_h5_compatible(metadata))
                    else:
                        hf[key].attrs.update(_h5_compatible(metadata))

        else:
            raise NotImplementedError(self.engine)

    def read_metadata(self, fname: str, key="default") -> dict:
        """
        Other Parameters
        ----------------
        key: str
            group where to read metadat from. If ``None``, add at root level. If ``'default'``,
            use engine's default (`/table` for `'vaex'`, `df` for `pytables`,
            root for `h5py` )
        """

        if self.engine in ["pytables", "pytables-fixed"]:
            fname = expanduser(fname)
            if key == "default":
                key = "df"
            with pd.HDFStore(fname, mode="r", complib="blosc", complevel=9) as f:
                try:
                    metadata = f.get_storer(key).attrs.metadata
                except KeyError as err:
                    print(f"Error reading metadata from {fname}")
                    raise err

        elif self.engine == "h5py":
            fname = expanduser(fname)
            if key == "default":
                key = None
            with h5py.File(fname, "r") as hf:
                if key is None:  # read metadta at root level
                    metadata = dict(hf.attrs)
                else:
                    try:
                        metadata = dict(hf[key].attrs)
                    except KeyError as err:
                        print(f"Error reading metadata from {fname}")
                        raise err

        elif self.engine == "vaex":
            if key == "default":
                key = r"/table"
            if isinstance(fname, list):
                metadata = []
                for f in fname:
                    f = expanduser(f)
                    with h5py.File(f, "r") as hf:
                        if key is None:  # read metadta at root level
                            metadata.append(dict(hf.attrs))
                        else:
                            metadata.append(dict(hf[key].attrs))
            else:
                fname = expanduser(fname)
                with h5py.File(fname, "r") as hf:
                    if key is None:  # add metadta at root level
                        metadata = dict(hf.attrs)
                    else:
                        try:
                            metadata = dict(hf[key].attrs)
                        except KeyError as err:
                            print(f"Error reading metadata from {fname}")
                            raise err

        else:
            raise NotImplementedError(
                f"'{self.engine}' is not implemented. Use 'pytables' or 'vaex' ?"
            )

        return metadata

    @classmethod
    def guess_engine(self, file, verbose=True):
        """Guess which HDF5 library ``file`` is compatible with

        .. note::
            it still take about 1 ms for this functino to execute. For extreme
            performance you want to directly give the correct engine

        Examples
        --------
        ::

            file = 'CO.hdf5'
            from radis.io.hdf5 import HDF5Manager
            engine = HDF5Manager.guess_engine(file)
            mgr = HDF5Manager(engine)
            mgr.read_metadata(file)

        """
        # See if it looks like PyTables
        import tables

        if tables.is_pytables_file(file):
            engine = "pytables"
        else:
            # Try Vaex
            file = expanduser(file)
            with h5py.File(file, mode="r") as hf:
                try:
                    hf[r"/table"]
                except KeyError:
                    engine = "h5py"
                else:
                    engine = "vaex"
        if verbose:
            print(f"Guessed that {file} was compatible with `{engine}` hdf5 engine")
        # raise
        return engine


def hdf2df(
    fname,
    columns=None,
    isotope=None,
    load_wavenum_min=None,
    load_wavenum_max=None,
    verbose=True,
    store_kwargs={},
    engine="guess",
):
    """Load a HDF5 line databank into a Pandas DataFrame.

    Adds HDF5 metadata in ``df.attrs``

    Parameters
    ----------
    fname : str
        HDF5 file name
    columns: list of str
        list of columns to load. If ``None``, returns all columns in the file.
    isotope: str
        load only certain isotopes : ``'2'``, ``'1,2'``, etc. If ``None``, loads
        everything. Default ``None``.
    load_wavenum_min, load_wavenum_max: float (cm-1)
        load only specific wavelength.

    Other Parameters
    ----------------
    store_kwargs: dict
        arguments forwarded to :py:meth:`~pandas.io.pytables.read_hdf`
    engine: ``'h5py'``, ``'pytables'``, ``'vaex'``, ``'auto'``
        which HDF5 library to use. If ``'guess'``, try to guess. Note: ``'vaex'``
        uses ``'h5py'`` compatible HDF5.

    Returns
    -------
    df: pandas Dataframe
        dataframe containing all lines or energy levels

    Examples
    --------

    ::


        path = getDatabankEntries("HITEMP-OH")['path'][0]
        df = hdf2df(path)

        df = hdf2df(path, columns=['wav', 'int'])

        df = hdf2df(path, isotope='2')
        df = hdf2df(path, isotope='1,2)

        df = hdf2df(path, load_wavenum_min=2300, load_wavenum_max=2500)

    Notes
    -----

    DataFrame metadata in ``df.attrs`` is still experimental in Pandas and can be lost
    during ``groupby, pivot, join or loc`` operations on the Dataframe.
    See https://stackoverflow.com/questions/14688306/adding-meta-information-metadata-to-pandas-dataframe

    Always check for existence !

    """

    if engine == "guess":
        engine = HDF5Manager.guess_engine(fname)

    t0 = time()

    # Selection
    selection = False
    if engine == "pytables":
        # Selection
        selection = True
        where = []
        if load_wavenum_min is not None:
            where.append(f"wav > {load_wavenum_min}")
        if load_wavenum_max is not None:
            where.append(f"wav < {load_wavenum_max}")
        if isotope:
            where.append(f'iso in {isotope.split(",")}')

    elif engine == "vaex":
        # Selection is done after opening the file time in vaex
        # see end of this function
        where = None
    else:
        raise NotImplementedError(engine)

    # Load :
    manager = HDF5Manager(engine)
    df = manager.load(fname, columns=columns, where=where, **store_kwargs)

    #  Selection in vaex
    if engine == "vaex":

        # Selection
        selection = True
        b = True
        if load_wavenum_min is not None:
            b *= df.wav > load_wavenum_min
        if load_wavenum_max is not None:
            b *= df.wav < load_wavenum_max
        if isotope is not None:
            from radis.misc.basics import is_float

            if is_float(isotope):
                b *= df.iso == int(isotope)
            else:
                b2 = False
                for iso in isotope.split(","):
                    b2 += df.iso == int(iso)
                b *= b2
        if b != True and False in b:
            df = df[b]  # note that this is a vaex Expression, not the DataFrame yet

        # Load
        if columns:  # load only these columns (if they exist)
            columns = [c for c in columns if c in df.columns]
        df = df.to_pandas_df(column_names=columns)

    # Read and add metadata in the DataFrame
    metadata = manager.read_metadata(fname)

    # Sanity Checks if loading the full file
    if not selection:
        if "total_lines" in metadata:
            assert len(df) == metadata["total_lines"]
        if "wavenumber_min" in metadata:
            assert df["wav"].min() == metadata["wavenumber_min"]
        if "wavenumber_max" in metadata:
            assert df["wav"].max() == metadata["wavenumber_max"]

    if isinstance(metadata, list):
        metadata_dict = {}
        for k, v in metadata[0].items():
            metadata_dict[k] = [v] + [M[k] for M in metadata[1:]]
        metadata = metadata_dict
    df.attrs.update(metadata)

    if verbose >= 3:
        from radis.misc.printer import printg

        printg(
            f"Generated dataframe from {fname} in {time()-t0:.2f}s ({len(df)} rows, {len(df.columns)} columns, {sys.getsizeof(df)*1e-6:.2f} MB)"
        )

    return df


#%%

if __name__ == "__main__":

    import pytest

    print("Testing factory:", pytest.main(["../test/io/test_hdf5.py"]))
