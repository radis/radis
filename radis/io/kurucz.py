""""

Summary
-----------------------

Kurucz database parser

-----------------------

Defines :func:`~radis.io.fetch_kurucz` based on :class:`~radis.api.kuruczapi.KuruczDatabaseManager`

"""

from os.path import abspath, expanduser, join

import radis
from radis.api.kuruczapi import KuruczDatabaseManager
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
    clean_cache_files=True,
    return_local_path=False,
    engine="default",
    output="pandas",
    parallel=True,
):
    """
    See e.g. :func:`~radis.io.hitemp.fetch_hitemp` for an explanation of the parameters largely applicable to `fetch_kurucz`

    .. note::

        if a registered entry already exists and `radis.config["ALLOW_OVERWRITE"]` is `True`:

        - if any situation arises where the databank needs to be re-downloaded, the possible urls are attempted in their usual order of preference, as if the databank hadn't been registered, rather than directly re-downloading from the same url that was previously registered, in case e.g. a new linelist has been uploaded since the databank was previously registered

        - If no partition function file is registered, e.g because one wasn't available server-side when the databank was last registered, an attempt is still made again to download it, to account for e.g. the case where one has since been uploaded
    """

    # largely based on :py:func:`~radis.io.fetch_geisa`
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
        parallel=parallel,
        engine=engine,
    )

    local_files, urlnames = [], []
    pf_path = []
    if ldb.is_registered():
        entries = getDatabankEntries(ldb.name)
        local_files, urlnames = entries["path"], entries["download_url"]
        if "parfunc" in entries:
            pf_path = [entries["parfunc"]]

    get_pf_files = True
    if ldb.is_registered() and not radis.config["ALLOW_OVERWRITE"]:
        error = False
        if cache == "regen" or not local_files:
            error = True
        files_to_check = local_files + pf_path
        if ldb.get_missing_files(files_to_check):
            error = True
        if not pf_path:
            get_pf_files = False  # assume partition function file is either unavailable or undesired for current species
        if error:
            raise Exception(
                'Changes are required to the local database, and hence updating the registered entry, but "ALLOW_OVERWRITE" is False. Set `radis.config["ALLOW_OVERWRITE"]=True` to allow the changes to be made and config file to be automatically updated accordingly.'
            )

    # Delete files if needed:

    if cache == "regen":
        ldb.remove_local_files(local_files + pf_path)
    ldb.check_deprecated_files(
        ldb.get_existing_files(local_files),
        auto_remove=True if cache != "force" else False,
    )

    get_main_files = True

    if len(local_files) > 1 or len(urlnames) > 1:
        raise Exception("only 1 database file is expected")

    if local_files and not ldb.get_missing_files(local_files):
        get_main_files = False
        ldb.actual_file = local_files[0]  # for ldb.load below
    if pf_path and not ldb.get_missing_files(pf_path):
        get_pf_files = False
        ldb.pf_path = pf_path[0]

    # Download files
    if get_main_files:
        main_files, main_urls = ldb.get_possible_files()
        for i in range(len(main_urls)):
            url = main_urls[i]
            file = main_files[i]
            print(f"Attempting to download {url}")
            try:
                ldb.download_and_parse([url], [file], 1)
            except OSError:
                if i == len(main_urls) - 1:  # all possible urls exhausted
                    print(f"Error downloading {url}.")
                    print(f"No source found for {ldb.molecule}")
                    raise
                else:
                    print(f"Error downloading {url}")
                    continue
            else:
                print(f"Successfully downloaded {url}")
                ldb.actual_file = file
                ldb.actual_url = url
                # local_files = [file]
                break  # no need to search any further

    if get_pf_files:
        pf_path, pf_url = ldb.get_pf_path()
        ldb.pf_path = pf_path
        try:
            ldb.download_and_parse([pf_url], [pf_path], 1)
        except OSError:
            print("a partition function file specific to this species was not found")
            get_pf_files = False
            ldb.pf_path = None

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
