""""

Summary
-----------------------

Kurucz database parser

-----------------------

Defines :py:func:`~radis.io.fetch_kurucz` based on :py:class:`~AdBKurucz`

"""

# import os
from os.path import abspath, expanduser, join

import radis
# from radis.api.hdf5 import DataFileManager
from radis.api.kuruczapi import KuruczDatabaseManager #AdBKurucz,get_atomic_number,get_ionization_state
from radis.misc.config import getDatabankEntries

def fetch_kurucz(
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
    potential_lowering=None
):
    """
    Note: if a registered entry already exists and `radis.config["ALLOW_OVERWRITE"]` is True:
    - if any situation arises where the databank needs to be re-downloaded, the possible urls are attempted in their usual order of preference, as if the databank hadn't been registered, rather than directly re-downloading from the same url that was previously registered, in case e.g. a new linelist has been uploaded since the databank was previously registered
    - If no partition function file is registered, e.g because one wasn't available server-side when the databank was last registered, an attempt is still made again to download it, to account for e.g. the case where one has since been uploaded
    """

    #largely based on :py:func:`~radis.io.fetch_geisa`
    if r"{molecule}" in databank_name:
        databank_name = databank_name.format(**{"molecule": molecule})

    if local_databases is None:

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

    local_files, urlnames = [], []
    pf_path = []
    if ldb.is_registered():
        entries = getDatabankEntries(ldb.name)
        local_files, urlnames = entries['path'], entries['download_url']
        if 'parfunc' in entries:
            pf_path = [entries['parfunc']]

    # # Get list of all expected local files for this database:
    # try:
    #     local_files, urlnames = ldb.get_filenames(return_reg_urls=True) # only expecting 1 file per molecule
    # except NotImplementedError:
    #     # no file registered so not sure what to expect
    #     local_files = []
    #     urlnames = []

    # main_files, main_urls = ldb.get_possible_files()
    # pf_path, pf_url = ldb.get_pf_path()
    # url_file = dict(zip(urlnames, local_files))
    
    get_pf_files = True
    if ldb.is_registered() and not radis.config["ALLOW_OVERWRITE"]:
        error = False
        if cache == "regen" or not local_files:
            error = True
        # main_url = list(set(main_urls) & set(urlnames))
        # if not main_url:
        #     error = True
        files_to_check = local_files + pf_path
        # files_to_check = [url_file[main_url[0]]]
        # if pf_url in urlnames:
        #     files_to_check.append(url_file[pf_url])
        # else:
        #     get_pf_files = False #assume partition function file is either unavailable or undesired for current species
        if ldb.get_missing_files(files_to_check):
            error = True
        if not pf_path:
            get_pf_files = False #assume partition function file is either unavailable or undesired for current species
        if error:
            raise Exception('Changes are required to the local database, and hence updating the registered entry, but "ALLOW_OVERWRITE" is False')
    
    # Delete files if needed:

    if cache == "regen":
        ldb.remove_local_files(local_files + pf_path)
    ldb.check_deprecated_files(
        ldb.get_existing_files(local_files),
        auto_remove=True if cache != "force" else False,
    )

    #download_files = []
    
    # missing_files = ldb.get_missing_files(local_files)
    # file_url = dict(zip(local_files, urlnames))
    # missing_urls = [file_url[file] for file in missing_files]
    
    # missing_url_file = dict(zip(missing_urls, missing_files))
    
    # main_url_file = dict(zip(main_urls, main_files))
    
    get_main_files = True

    if len(local_files) > 1 or len(urlnames) > 1:
        raise Exception('only 1 database file is expected')
    
    if local_files and not ldb.get_missing_files(local_files):
        get_main_files = False
        ldb.actual_file = local_files[0] # for ldb.load below
    if pf_path and not ldb.get_missing_files(pf_path):
        get_pf_files = False
        ldb.pf_path = pf_path[0]
    
    # for url in urlnames:
    #     if url not in missing_urls:
    #         if url in main_urls:
    #             get_main_files = False
    #             ldb.actual_url = url # in case of re-registering ...
    #             ldb.actual_file = url_file[url] # ... and for ldb.load below
    #         elif url == pf_url:
    #             get_pf_files = False
    #             ldb.pf_url = url
    #             ldb.pf_path = url_file[url]



    # for url in main_urls:
    #     if url not in urlnames or url in missing_urls:
    #         get_main_files = True
    #     else:
    #         ldb.actual_url = url
    #         ldb.actual_file = main_url_file[url]
    # if pf_url not in urlnames or pf_url in missing_urls:
    #     get_pf = True
    # else
    # for url in main_urls + pf_url:
    #     if url not in urlnames or url in missing_urls:
    #         if url in main_urls:
    #             get_main_files = True
    #     elif url in main_urls: #in case it gets (re-)registered, then just register at this url and path
    #         ldb.actual_url = url
    #         ldb.actual_file = main_url_file[url]
        
    
    
    # if missing_files:
    #     file_url = dict(zip(local_files, urlnames))
    #     missing_urls = [file_url[file] for file in missing_files]
    #     main_files, main_urls = ldb.get_possible_files()
    #     pf_path, pf_url = ldb.get_pf_path()
    #     main_url_file = dict(zip(main_urls, main_files))
    #     missing_url_file = dict(zip(missing_urls, missing_files))
    #     for url in missing_urls:
    #         if url in main_urls:
    #             missing_url_file.pop(url)
    #             get_main_files = True
    #         if url == pf_url:
    #             missing_url_file[url] = pf_path
        
    
    # existing_files = ldb.get_existing_files(main_files)
    # if len(existing_files) == 1:
    #     ldb.actual_file = existing_files[0]
    #     main_file_url = dict(zip(main_files, main_urls))
    #     download_files = []
    # # elif urlnames is not None:
    # #     if not local_files:
    # #         download_files, _ = ldb.get_possible_files(urlnames=urlnames)
    # #     else:
    # #         download_files = local_files
    # else:
    #     download_files, download_urls = ldb.get_possible_files()
    
    # Download files
    if get_main_files:
        main_files, main_urls = ldb.get_possible_files()
        # if urlnames is None:
        #     urlnames = ldb.fetch_urlnames()
        # filesmap = dict(zip(download_files, urlnames))
        # download_urls = [filesmap[k] for k in download_files]
        for i in range(len(main_urls)):
            url = main_urls[i]
            file = main_files[i]
            print(f'Attempting to download {url}')
            try:
                ldb.download_and_parse([url], [file], 1)
            except OSError:
                if i == len(main_urls) - 1: #all possible urls exhausted
                    print(f'Error downloading {url}.')
                    print(f'No source found for {ldb.molecule}')
                    raise
                else:
                    print(f'Error downloading {url}')
                    continue
            else:
                print(f'Successfully downloaded {url}')
                ldb.actual_file = file
                ldb.actual_url = url
                # local_files = [file]
                break #no need to search any further
    
    if get_pf_files:
        pf_path, pf_url = ldb.get_pf_path()
        ldb.pf_path = pf_path
        # ldb.pf_url = pf_url
        try:
            ldb.download_and_parse([pf_url], [pf_path], 1)
        except OSError:
            print('a partition function file specific to this species was not found')
            get_pf_files = False
            ldb.pf_path = None
            # ldb.pf_url = None

    # Register
    if get_main_files or get_pf_files or not ldb.is_registered():
        ldb.register(get_main_files, get_pf_files)

    if (get_main_files or get_pf_files) and clean_cache_files:
        ldb.clean_download_files()

    if isotope and type(isotope) == int:
        isotope = str(isotope)

    # Load and return
    df = ldb.load(
        [ldb.actual_file],
        columns=columns,
        within=[("iso", isotope)] if isotope is not None else [],
        # for relevant files, get only the right range :
        lower_bound=[("wav", load_wavenum_min)] if load_wavenum_min is not None else [],
        upper_bound=[("wav", load_wavenum_max)] if load_wavenum_max is not None else [],
        output=output,
    )

    # based on ExoMol:
    if output in ["pandas", "vaex"]:  # no attribtes in "Jax" or "Vaex" mode
        attrs = {}
        attrs["molecule"] = molecule

        if output == "vaex":
            df.attrs = {}
            df.attrs = attrs
        elif output == "pandas":
            for k, v in attrs.items():
                df.attrs[k] = v

    return (df, local_files, ldb.pf_path) if return_local_path else df



# def fetch_kurucz(species, isotope, engine):
#     kurucz = AdBKurucz(species)
#     data_manager = DataFileManager(engine=engine)
#     atomic_number = f"{get_atomic_number(species):02}"
#     ionization_state_str = f"{get_ionization_state(species):02}"
#     hdf5_file = f"gf{atomic_number}{ionization_state_str}.hdf5"
#     kurucz.url = kurucz.get_url(atomic_number, ionization_state_str)
#     # kurucz.hdf5_file = hdf5_file  # Set kurucz.hdf5_file to hdf5_file

#     # If hdf5 file exists, read data from it
#     if os.path.exists(hdf5_file):
#         print("HDF5 file already exists, reading data from it.")
#     else:
#         kuruczf = kurucz.download_file()
#         df = kurucz.read_kurucz(kuruczf)
#         data_manager.write(
#             file=hdf5_file,
#             df=df,
#             append=False,
#             format="table",
#             data_columns=df.columns,
#         )
#     df = data_manager.load(fname=hdf5_file, within=[("iso", isotope)] if isotope is not None else [])
#     #kurucz.add_airbrd(df)

#     return hdf5_file, df
