""""

Summary
-----------------------

Kurucz database parser

-----------------------

Defines :py:func:`~radis.io.fetch_kurucz` based on :py:class:`~AdBKurucz`

"""

import os
from os.path import abspath, expanduser, join

from radis.api.hdf5 import DataFileManager
from radis.api.kuruczapi import AdBKurucz,get_atomic_number,get_ionization_state, KuruczDatabaseManager

def fetch_kurucz_new(
    molecule,
    local_databases=None,
    databank_name="Kurucz-{molecule}",
    isotope=None,
    load_wavenum_min=None,
    load_wavenum_max=None,
    columns=None,
    cache=True,
    verbose=True,
    chunksize=100000,
    clean_cache_files=True,
    return_local_path=False,
    engine="default",
    output="pandas",
    parallel=True,
):
    #largely based on :py:func:`~radis.io.fetch_geisa`
    if r"{molecule}" in databank_name:
        databank_name = databank_name.format(**{"molecule": molecule})

    if local_databases is None:
        import radis

        local_databases = join(radis.config["DEFAULT_DOWNLOAD_PATH"], "kurucz")
    local_databases = abspath(expanduser(local_databases))

    ldb = KuruczDatabaseManager(
        databank_name,
        molecule=molecule,
        local_databases=local_databases,
        verbose=verbose,
        chunksize=chunksize,
        parallel=parallel,
        engine=engine,
    )

    # Get list of all expected local files for this database:
    try:
        local_files, _ = ldb.get_filenames(return_reg_urls=True) # only expecting 1 file per molecule
    except NotImplementedError:
        # no file registered so not sure what to expect
        local_files = []
        # urlnames = None

    # Delete files if needed:

    if cache == "regen":
        ldb.remove_local_files(local_files)
    ldb.check_deprecated_files(
        ldb.get_existing_files(local_files),
        auto_remove=True if cache != "force" else False,
    )

    existing_files = ldb.get_existing_files(local_files)
    if existing_files:
        download_files = []
    # elif urlnames is not None:
    #     if not local_files:
    #         download_files, _ = ldb.get_possible_files(urlnames=urlnames)
    #     else:
    #         download_files = local_files
    else:
        download_files, download_urls = ldb.get_possible_files()
    
    # Download files
    if len(download_files) > 0:
        # if urlnames is None:
        #     urlnames = ldb.fetch_urlnames()
        # filesmap = dict(zip(download_files, urlnames))
        # download_urls = [filesmap[k] for k in download_files]
        for i in range(len(download_urls)):
            url = download_urls[i]
            file = download_files[i]
            print(f'Attempting to download {url}')
            try:
                ldb.download_and_parse([url], [file], 1)
            except OSError:
                if i == len(download_urls) - 1:
                    print(f'Error downloading {url}.')
                    print(f'No source found for {ldb.molecule}')
                    raise
                else:
                    print(f'Error downloading {url}. Attempting to download {download_urls[i+1]}')
                    continue
            else:
                print(f'Successfully downloaded {url}')
                ldb.actual_file = file
                ldb.actual_url = url
                local_files = [file]
                break

    # Register
    if not ldb.is_registered():
        ldb.register()

    if len(download_files) > 0 and clean_cache_files:
        ldb.clean_download_files()

    if isotope and type(isotope) == int:
        isotope = str(isotope)

    # Load and return
    df = ldb.load(
        local_files,  # filter other files,
        columns=columns,
        within=[("iso", isotope)] if isotope is not None else [],
        # for relevant files, get only the right range :
        lower_bound=[("wav", load_wavenum_min)] if load_wavenum_min is not None else [],
        upper_bound=[("wav", load_wavenum_max)] if load_wavenum_max is not None else [],
        output=output,
    )

    return (df, local_files) if return_local_path else df



def fetch_kurucz(species, isotope, engine):
    kurucz = AdBKurucz(species)
    data_manager = DataFileManager(engine=engine)
    atomic_number = f"{get_atomic_number(species):02}"
    ionization_state_str = f"{get_ionization_state(species):02}"
    hdf5_file = f"gf{atomic_number}{ionization_state_str}.hdf5"
    kurucz.url = kurucz.get_url(atomic_number, ionization_state_str)
    # kurucz.hdf5_file = hdf5_file  # Set kurucz.hdf5_file to hdf5_file

    # If hdf5 file exists, read data from it
    if os.path.exists(hdf5_file):
        print("HDF5 file already exists, reading data from it.")
    else:
        kuruczf = kurucz.download_file()
        df = kurucz.read_kurucz(kuruczf)
        data_manager.write(
            file=hdf5_file,
            df=df,
            append=False,
            format="table",
            data_columns=df.columns,
        )
    df = data_manager.load(fname=hdf5_file, within=[("iso", isotope)] if isotope is not None else [])
    #kurucz.add_airbrd(df)

    return hdf5_file, df
