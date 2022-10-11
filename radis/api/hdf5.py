# -*- coding: utf-8 -*-
"""
Defines the :py:class:`~radis.api.hdf5.DataFileManager` class
"""

import os
import pathlib
import sys
from os.path import abspath, exists, expanduser, splitext
from time import time

import h5py
import pandas as pd
from tables.exceptions import NoSuchNodeError


def vaexsafe_colname(name):
    """replace '/' (forbidden in HDF5 vaex column names with '_'
    https://github.com/radis/radis/issues/473
    https://github.com/vaexio/vaex/issues/1255
    """
    return name.replace("/", "_")


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

    pytables_manager = DataFileManager(engine="pytables")

    # Read metadata
    file_metadata = pytables_manager.read_metadata(fname)

    # Write Vaex file
    df.export_hdf5(fname_vaex)
    df.close()  # try to fix file not closed()  TODO: remove?
    del df  # same TODO

    vaex_manager = DataFileManager(engine="vaex")
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
    def __init__(*args, **kwargs):
        raise DeprecationWarning("HDF5Manager replaced with DataFileManager")


class DataFileManager(object):
    def __init__(self, engine=None):
        """Class to handle all DataFrame-librairies with one common API

        All functions may not be fully implemetned, will raise a NotImplementedError
        if that's not the case.

        Librairies ::

            'vaex'     > HDF5,  column-based
            'pytables' > Pandas's HDF5,  row-based
            'h5py'     > HDF5
            'feather'  > feather

        Functions ::

            add_metadata
            read_metadata
            write
            read
            guess_engine

        Examples
        --------
        ::

            file = 'CO.hdf5'
            from radis.api.hdf5 import DataFileManager
            engine = DataFileManager.guess_engine(file)
            mgr = DataFileManager(engine)
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
        append=False,
        key="default",
        format="table",
        data_columns=["iso", "wav", "nu_lines"],
    ):
        """Write dataframe ``df`` to ``file``

        Parameters
        ----------
        df: DataFrame

        Other Parameters
        ----------------
        key: str
            group to write to. If ``None``, write at root level. If ``'default'``,
            use engine's default (`/table` for `'vaex'`, `df` for `pytables`,
            root for `h5py` )
        data_columns : list
            only these column names will be searchable directly on disk to
            load certain lines only. See :py:func:`~radis.api.hdf5.hdf2df`
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
            if isinstance(df, pd.DataFrame):
                import vaex

                df = vaex.from_pandas(df)

            for c in df.columns:  # remove "/" in columns (forbidden)
                df.rename(c, vaexsafe_colname(c))

            if key == "default":
                key = r"/table"

            if append == True:
                # In vaex we cannot append. Here we write lots of small files then combine them.
                # self.combine_temp_batch_files() should be called at the end.

                # To start with, get an available temp file name :
                base, ext = splitext(file)
                i = 0
                temp_batch_file = base + "_temp" + str(i).zfill(5) + ext
                while temp_batch_file in self._temp_batch_files:
                    i += 1
                    temp_batch_file = base + "_temp" + str(i).zfill(5) + ext
                file = temp_batch_file
                # Check no remaining one from a non-cleaned previous run:
                if exists(file):
                    from radis.misc.printer import printr

                    printr(f"Temp file {temp_batch_file} already exists: deleting it")
                    os.remove(file)
                self._temp_batch_files.append(file)
            # Write:
            df.export_hdf5(file, group=key, mode="w")
        elif self.engine == "feather":
            df.to_feather(file)
        else:
            raise NotImplementedError(self.engine)
            # h5py is not designed to write Pandas DataFrames

    def get_columns(self, local_file):
        """Get all columns (without loading all Dataframe)"""
        engine = self.engine
        local_file = expanduser(local_file)
        if engine == "vaex":
            import vaex

            # by default vaex does not load everything
            df = vaex.open(local_file)
            columns = df.columns
            df.close()

        elif engine == "pytables":
            with pd.HDFStore(local_file, "r") as store:
                columns = store.select("df", start=1, stop=1).columns
        elif engine in ["h5py"]:
            raise NotImplementedError
        else:
            raise ValueError(engine)

        return columns

    def combine_temp_batch_files(
        self, file, key="default", sort_values=None, delete_nan_columns=True
    ):
        """Combine all batch files in ``self._temp_batch_files`` into one.
        Removes all batch files.
        """
        file = expanduser(file)
        if self.engine == "vaex":
            if len(self._temp_batch_files) == 0:
                # No temp file created. File is probably already created (append=False mode)
                # Else, something unexpected happens --> raise error
                if exists(file):
                    return
                raise ValueError(f"No batch temp files were written for {file}")
            if key == "default":
                key = r"/table"
            import vaex

            df = vaex.open(self._temp_batch_files, group=key)
            # Removing Nan values columns
            if delete_nan_columns:
                import numpy as np

                for column in df.columns:
                    col = df[column].values
                    if type(col[0]) in [np.int32, np.float64] and np.isnan(np.sum(col)):
                        del df[column]

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
        """clean before deleting"""
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
        lower_bound=[],
        upper_bound=[],
        within=[],
        output="pandas",
        **store_kwargs,
    ):
        """
        Other Parameters
        ----------------
        columns: list of str
            list of columns to load. If ``None``, returns all columns in the file.
        output: 'pandas', 'vaex', 'jax'
            format of the output DataFrame. If ``'jax'``, returns a dictionary of
            jax arrays.
        lower_bound: list of tuples [(column, lower_bound), etc.]
            ::

                lower_bound =[("wav", load_wavenum_min)]
        upper_bound_bound: list of tuples [(column, upper_bound), etc.]
            ::

                upper_bound=[("wav", load_wavenum_max)]
        within: list of tuples [(column, within_list), etc.]
            ::

                within=[("iso", isotope.split(","))]
        """
        if not lower_bound and not upper_bound and not within:
            df = self.read(
                fname,
                columns,
                **store_kwargs,
            )
        else:
            df = self.read_filter(
                fname,
                columns,
                lower_bound=lower_bound,
                upper_bound=upper_bound,
                within=within,
                **store_kwargs,
            )

        # Load
        # Convert from reading format ("engine") to output format ("output")
        engine = self.engine

        if engine == output:
            df = df
        elif engine == "vaex":
            # in vaex, column selection has to happen now
            if columns:  # load only these columns (if they exist)
                columns = [c for c in columns if c in df.columns]
            if output == "pandas":
                df_pandas = df.to_pandas_df(column_names=columns)
                df.close()
                df = df_pandas
            elif output == "jax":
                if columns == None:
                    columns = list(df.columns)
                out = {}
                try:
                    import jax.numpy as jnp
                except ImportError:
                    print("Jax not found. Using Numpy.")
                    import numpy as jnp
                for c in columns:
                    if (
                        c == "Sij0"
                    ):  # (special case for Sij0 : we store logSij0 instead, to store them as jax np.32 arrays)
                        import numpy as np

                        out["logsij0"] = jnp.array(np.log(df[c].values))
                    else:
                        out[c] = jnp.array(df[c].values)
                df.close()
                df = out
            else:
                raise NotImplementedError(f"output {output} for engine {engine}")
        elif engine == "pytables":
            if output == "pandas":
                df = df
            else:
                raise NotImplementedError(f"output {output} for engine {engine}")
        elif engine == "feather":
            if output == "pandas":
                df = df
            else:
                raise NotImplementedError(f"output {output} for engine {engine}")
        else:
            raise NotImplementedError(output)

        if output == "vaex":
            df = df.extract()  # return DataFrame containing only the filtered rows

        return df

    def read(
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
        fname: str
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
                        if f.endswith(".hdf5") and exists(f.replace(".hdf5", ".h5")):
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

        elif self.engine == "feather":
            assert where is None
            fname = expanduser(fname)
            return pd.read_feather(fname)

        else:
            raise NotImplementedError(self.engine)

        return df

    def read_filter(
        self,
        fname,
        columns=None,
        lower_bound=[],
        upper_bound=[],
        within=[],
        **store_kwargs,
    ):
        """
        Parameters
        ----------
        fname: str
        columns: list of str
            list of columns to load. If ``None``, returns all columns in the file.
        lower_bound: list of tuples [(column, lower_bound), etc.]
            ::

                lower_bound =[("wav", load_wavenum_min)]
        upper_bound_bound: list of tuples [(column, upper_bound), etc.]
            ::

                upper_bound=[("wav", load_wavenum_max)]
        within: list of tuples [(column, within_list), etc.]
            ::

                within=[("iso", isotope.split(","))]

        """

        # Selection
        if self.engine == "pytables":
            # Selection
            where = []
            for (column, lbound) in lower_bound:
                where.append(f"{column} > {lbound}")
            for (column, ubound) in upper_bound:
                where.append(f"{column} < {ubound}")
            for (column, withinv) in within:
                where.append(f"{column} in {withinv.split(',')}")

        elif self.engine in ["vaex", "feather"]:
            # Selection is done after opening the file time in vaex
            # see end of this function
            where = None
        else:
            raise NotImplementedError(self.engine)

        # Load :
        df = self.read(fname, columns=columns, where=where, **store_kwargs)

        #  Selection in vaex
        if self.engine in ["vaex", "feather"]:
            # (note that in Vaex, the selection happens on disk whereas Feather
            # is already loaded as a Pandas DataFrame in RAM)

            # Selection
            b = True
            for (column, lbound) in lower_bound:
                b *= df[column] > lbound
            for (column, ubound) in upper_bound:
                b *= df[column] < ubound
            for (column, withinv) in within:
                b2 = False
                for val in withinv.split(","):
                    b2 += df[column] == float(val)
                b *= b2
            if b is not True and False in b:
                df = df[
                    b
                ]  # note in Vaex mode, this is a vaex Expression, not the DataFrame yet

        return df

    def cache_file(self, fname):
        """Return the corresponding cache file name for fname.

        Other Parameters
        ----------------
        engine: ``'h5py'``, ``'pytables'``, ``'vaex'``
           which HDF5 library to use. Default ``pytables``
        """
        if self.engine in ["pytables", "pytables-fixed"]:
            return pathlib.Path(fname).with_suffix(".h5")
        elif self.engine in ["h5py", "vaex"]:
            return pathlib.Path(fname).with_suffix(".hdf5")
        elif self.engine == "feather":
            return pathlib.Path(fname).with_suffix(".feather")
        else:
            raise ValueError(self.engine)

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
        from radis.api.cache_files import _h5_compatible

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
                except AttributeError as err:
                    if "Attribute 'metadata' does not exist" in str(err):
                        metadata = {}
                    else:
                        raise err

        elif self.engine == "feather":
            return {}  # no metadata

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
                        except (KeyError, OSError) as err:
                            if key == r"/table":
                                # backward compat : some old files were generated with key=None
                                metadata = dict(hf.attrs)
                                print(
                                    f"Error reading metadata from {fname}. Regenerate file one day?"
                                )
                            else:
                                print(f"Error reading metadata from {fname}")
                                raise err

        else:
            raise NotImplementedError(
                f"'{self.engine}' is not implemented. Use 'pytables' or 'vaex' ?"
            )

        return metadata

    def to_numpy(self, df):
        """Convert DataFrame to numpy"""

        if self.engine == "vaex":
            import vaex

            return vaex.array_types.to_numpy(df)
        elif self.engine == "feather":
            return df.to_numpy()
        else:
            raise NotImplementedError(self.engine)

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
        if file.endswith(".feather"):
            engine = "feather"
        else:
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
            print(f"Guessed that {file} was compatible with `{engine}` engine")
        # raise
        return engine

    def has_nan(self, column):
        if self.engine == "vaex":
            b = column.isnan()  # TODO: check if can be made faster?
            return b.sum() > 0
        elif self.engine in ["pytables", "feather"]:
            return column.hasnans
        else:
            raise NotImplementedError(self.engine)


def hdf2df(
    fname,
    columns=None,
    isotope=None,
    load_wavenum_min=None,
    load_wavenum_max=None,
    verbose=True,
    store_kwargs={},
    engine="guess",
    output="pandas",
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
    output: 'pandas', 'vaex', 'jax'
        format of the output DataFrame. If ``'jax'``, returns a dictionary of
        jax arrays.

    Returns
    -------
    df: pandas Dataframe, or vaex DataFrameLocal, or dictionary of Jax arrays
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

    fname = abspath(expanduser(fname))

    if engine == "guess":
        engine = DataFileManager.guess_engine(fname)

    t0 = time()

    manager = DataFileManager(engine)

    df = manager.load(
        fname,
        columns,
        lower_bound=[("wav", load_wavenum_min)] if load_wavenum_min is not None else [],
        upper_bound=[("wav", load_wavenum_max)] if load_wavenum_max is not None else [],
        within=[("iso", isotope)] if isotope is not None else [],
        output=output,
    )

    # Read and add metadata in the DataFrame
    metadata = manager.read_metadata(fname)

    # Sanity Checks if loading the full file
    selection = isotope or load_wavenum_min or load_wavenum_max
    if not selection:
        if "total_lines" in metadata:
            assert len(df) == metadata["total_lines"]
        if "wavenumber_min" in metadata:
            assert df["wav"].min() == metadata["wavenumber_min"]
        if "wavenumber_max" in metadata:
            assert df["wav"].max() == metadata["wavenumber_max"]

    # Update metadata
    if output != "jax":  # metadata not given in Jax output
        if isinstance(metadata, list):
            metadata_dict = {}
            for k, v in metadata[0].items():
                metadata_dict[k] = [v] + [M[k] for M in metadata[1:]]
            metadata = metadata_dict
        if (
            output in "vaex"
        ):  # by default vaex.dataframe.DataFrameLocal doesn't have a .attrs attribute
            df.attrs = {}
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
