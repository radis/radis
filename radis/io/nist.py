""" "

Summary
-----------------------

NIST database parser

-----------------------

Defines :func:`~radis.io.fetch_nist` based on :class:`~radis.api.nistapi.NISTDatabaseManager`

"""

from os.path import abspath, expanduser, join

import radis
from radis.api.nistapi import NISTDatabaseManager


def fetch_nist(
    molecule,  # replace
    local_databases=None,
    databank_name="NIST-{molecule}",
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
    See e.g. :func:`~radis.io.hitemp.fetch_hitemp` for an explanation of the parameters largely applicable to `fetch_nist`

    .. note::

        if a registered entry already exists and `radis.config["ALLOW_OVERWRITE"]` is `True`:

        - if any situation arises where the databank needs to be re-downloaded, the possible urls are attempted in their usual order of preference, as if the databank hadn't been registered, rather than directly re-downloading from the same url that was previously registered, in case e.g. a new linelist has been uploaded since the databank was previously registered

    """

    # largely based on :py:func:`~radis.io.fetch_geisa`
    if r"{molecule}" in databank_name:
        databank_name = databank_name.format(**{"molecule": molecule})

    if local_databases is None:

        local_databases = join(radis.config["DEFAULT_DOWNLOAD_PATH"], "NIST")
    local_databases = abspath(expanduser(local_databases))

    ldb = NISTDatabaseManager(
        databank_name,
        molecule=molecule,
        local_databases=local_databases,
        verbose=verbose,
        parallel=parallel,
        engine=engine,
    )

    # Get list of all expected local files for this database:
    local_files, urlnames = ldb.get_filenames()

    if ldb.is_registered() and not radis.config["ALLOW_OVERWRITE"]:
        error = False
        if cache == "regen" or not local_files:
            error = True
        if ldb.get_missing_files(local_files):
            error = True
        if error:
            raise Exception(
                'Changes are required to the local database, and hence updating the registered entry, but "ALLOW_OVERWRITE" is False. Set `radis.config["ALLOW_OVERWRITE"]=True` to allow the changes to be made and config file to be automatically updated accordingly.'
            )

    # Delete files if needed:

    if cache == "regen":
        ldb.remove_local_files(local_files)
    ldb.check_deprecated_files(
        ldb.get_existing_files(local_files),
        auto_remove=True if cache != "force" else False,
    )

    get_main_files = True

    if len(local_files) > 1 or len(urlnames) > 1:
        raise Exception("only 1 database file is expected")

    if local_files and not ldb.get_missing_files(local_files):
        get_main_files = False

    # Download files
    if get_main_files:
        urlnames = ldb.fetch_urlnames()
        ldb.download_and_parse(urlnames, local_files, 1)

    # Register
    if get_main_files or not ldb.is_registered():
        ldb.register(local_files, urlnames)

    if get_main_files and clean_cache_files:
        ldb.clean_download_files()

    if isotope and type(isotope) == int:
        isotope = str(isotope)

    # Load and return
    df = ldb.load(
        local_files,
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

    return (df, local_files) if return_local_path else df
