# -*- coding: utf-8 -*-
"""Implements a spectrum database :class:`~radis.tools.database.SpecDatabase`
class to manage them all.

It basically manages a list of Spectrum JSON files, adding a Pandas
dataframe structure on top to serve as an efficient index to visualize
the spectra input conditions, and slice through the Dataframe with
easy queries

Examples
--------

See and get objects from database::

    from radis.tools import SpecDatabase
    db = SpecDatabase(r"path/to/database")     # create or loads database

    db.update()  # in case something changed (like a file was added manually)
    db.see(['Tvib', 'Trot'])   # nice print in console

    s = db.get('Tvib==3000 & Trot==1500')[0]  # get all spectra that fit conditions
    db.add(s)  # update database (and raise error because duplicate!)

Note that :py:class:`~radis.lbl.factory.SpectrumFactory` can be configured to
automatically look-up and update a database when spectra are calculated.

An example of script to update all spectra conditions in a database (ex: when
a condition was added afterwards to the Spectrum class)::

    # Example: add the 'medium' key in conditions
    db = "database_CO"
    for f in os.listdir(db):
       if not f.endswith('.spec'): continue
       s = load_spec(join(db, f))
       s.conditions['medium'] = 'vacuum'
       s.store(join(db,f), if_exists_then='replace')

You can see more examples on the :ref:`Spectrum Database section <label_spectrum_database>`
of the website.

-------------------------------------------------------------------------------
"""

# TODO:

# - Alert if case already in database when generating from a SpectrumFactory
# connected to a SpecDatabase

import os
import sys
from os.path import (
    abspath,
    basename,
    dirname,
    exists,
    getsize,
    isdir,
    join,
    split,
    splitext,
)
from shutil import copy2
from time import strftime
from warnings import warn

import json_tricks
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from numpy import array
from scipy.interpolate import griddata

from radis.misc.basics import all_in, is_float, list_if_float
from radis.misc.debug import printdbg
from radis.misc.printer import printr
from radis.spectrum.spectrum import Spectrum

_scalable_inputs = ["mole_fraction", "path_length"]

# %% Tools


def is_jsonable(x):
    try:
        #        json.dumps(x)
        json_tricks.dumps(x)
        return True
    except:
        return False


# def jsonize(x):
#    ''' Converts x so it can be stored in JSON
#    Replaced by call to json-tricks '''
#    typ = type(x)
#
#    # Convert object based on type
#    if typ == dict:
#        # If dictionary, deal with all keys recursively
#        out = {}
#        for k, v in x.items():
#            out[k] = jsonize(v)
#    elif isinstance(typ, str):
#        # If it looks like string, store as raw text (that fixes most trouble with paths)
#        out = r'{0}'.format(x)
#    elif typ == pd.DataFrame:
#        out = x.to_json()
#    elif typ == np.ndarray:
#        # converts numpy array to list
#        out = x.tolist()
#    elif typ in [np.float64, np.float32]:
#        # cast numpy floats to float
#        out = float(x)
#    elif typ in [np.int32, np.int64]:
#        # cast numpy int to int
#        out = int(x)
#    else:
#        # Keep the object as it is
#        out = x
#
#    # Check it's jsonable
#    if is_jsonable(out):
#        return out
#    else:
#        if out is x:
#            raise TypeError('{0} (type {1}) is not jsonable'.format(out, type(out)))
#        else:
#            raise TypeError('{0} (type {1}) converted from {2} (type {3}) is not jsonable'.format(
#                    out, type(out), x, typ))

# %% Serializing functions
# ... functions to store / load a Spectrum to / from a JSON file


def save(
    s: Spectrum,
    path,
    discard=[],
    compress=True,
    add_info=None,
    add_date=None,
    if_exists_then="increment",
    verbose=True,
    warnings=True,
):
    """Save a :class:`~radis.spectrum.spectrum.Spectrum` object in JSON format.
    Object can be recovered with :func:`~radis.tools.database.load_spec`. If
    many :class:`~radis.spectrum.spectrum.Spectrum` are saved in a same folder
    you can view their properties with the
    :class:`~radis.tools.database.SpecDatabase` structure.

    Parameters
    ----------
    s: Spectrum
        to save
    path: str
        filename to save. No extension needed. If filename already
        exists then a digit is added. If filename is a directory then a new
        file is created within this directory.
    discard: list of str
        parameters to discard. To save some memory.
    compress: boolean
        if ``False``, save under text format, readable with any editor.
        if ``True``, saves under binary format. Faster and takes less space.
        If ``2``, removes all quantities that can be regenerated with s.update(),
        e.g, transmittance if abscoeff and path length are given, radiance if
        emisscoeff and abscoeff are given in non-optically thin case, etc.
        Default ``False``
    add_info: list, or None/False
        append these parameters and their values if they are in conditions.
        e.g::

            add_info = ['Tvib', 'Trot']

    add_date: str, or ``None``/``False``
        adds date in strftime format to the beginning of the filename.
        e.g::

            add_date = '%Y%m%d'

    if_exists_then: ``'increment'``, ``'replace'``, ``'error'``, ``'ignore'``
        what to do if file already exists. If ``'increment'`` an incremental digit
        is added. If ``'replace'`` file is replaced (!). If ``'ignore'`` the
        Spectrum is not added to the database and no file is created.
        If ``'error'`` (or anything else) an error is raised. Default ``'increment'``.

    Returns
    -------
    fout: str
        filename used (may be different from given path as new info or
        incremental identifiers are added)


    See Also
    --------
    :class:`~radis.tools.database.SpecDatabase`,
    :func:`~radis.tools.database.load_spec`,
    :meth:`~radis.spectrum.spectrum.Spectrum.store`
    """

    # 1) Format to JSON writable dictionary
    sjson = _format_to_jsondict(s, discard, compress, verbose=verbose)

    # 2) Get final output name (add info, extension, increment number if needed)
    fout = _get_fout_name(path, if_exists_then, add_date, add_info, sjson, verbose)
    if exists(fout) and if_exists_then == "ignore":
        return fout

    # 3) Now is time to save
    if compress:
        with open(fout, "wb") as f:
            json_tricks.dump(
                sjson,
                f,
                compression=True,  # calls gzip compression
                properties={
                    "ndarray_compact": True
                },  # use compact numpy format in json-tricks 3.15+
            )
    else:
        with open(fout, "w") as f:
            #        json.dump(sjson, f)
            json_tricks.dump(sjson, f, indent=4)  # takes more space but more readable

    if verbose:
        print(
            "Spectrum stored in {0} ({1:.1f}Mb)".format(
                fout, getsize(fout) / 1024 * 1e-3
            )
        )

    return fout  # return final name


def _format_to_jsondict(s: Spectrum, discard, compress, verbose=True):
    """Format to JSON writable dictionary.

    Notes
    -----

    path names create much troubles on reload if they are stored with '/'
    Make sure we use raw format. This is now dealt with by the json-tricks module

    We also store the initial type of all conditions so they are correctly
    reproduced on loading .
    """
    # TODO: convert functions to text with marshal:
    # see: https://stackoverflow.com/a/51938459/5622825

    # Add all in a dictionary using json-tricks (that deals with numpy and pandas
    # array, as well as text: see 24/03/18 for former implementation that was
    # manually converted to jsonable)
    sjson = {}
    for attr in s.__slots__:
        sjson[attr] = s.__getattribute__(attr)

    # Discard some entries
    for k in discard:
        if k in sjson:
            del sjson[k]

    # Check that all conditions are jsonable, discard if not
    for k in list(sjson["conditions"].keys()):
        v = sjson["conditions"][k]
        if not is_jsonable(v):
            if verbose:
                printr(
                    "Discarded {0} from conditions as not jsonable ({1})".format(
                        k, type(v)
                    )
                )
            del sjson["conditions"][k]

    # if compress>=2, remove unecessary spectral quantities (that can be recomputed
    # from the rest)
    if compress >= 2:
        sjson["_q"] = sjson["_q"].copy()
        sjson = _compress(s, sjson)

    return sjson


def _get_fout_name(path, if_exists_then, add_date, add_info, sjson, verbose):
    """Get final output name   (add info, extension, increment number if
    needed)"""

    conditions = sjson["conditions"]

    if isdir(path):
        fold, name = path, ""
    else:
        fold, name = split(path)

    # ... add date info
    if add_date not in ["", None, False]:
        date = strftime(add_date)
    else:
        date = ""

    # ... add conditions info
    if isinstance(add_info, str):
        add_info = [add_info]
    if add_info not in [[], {}, None, False]:
        # complete name with info about calculation conditions
        info = []
        for k in add_info:
            if k in conditions:
                v = conditions[k]
                # Format info
                # ... special cases
                if k in ["Tvib"] and is_float(v):
                    # (dev): it can happen in multi-Tvib mode that Tvib is
                    # a str ('Tvib1,Tvib2,Tvib3') and not a float
                    vs = "{0:.0f}".format(v)
                elif k in ["Telec", "Tgas", "Trot"]:
                    vs = "{0:.0f}".format(v)
                # ... general case
                elif is_float(v):
                    vs = "{0:.3g}".format(v)
                else:
                    vs = "{0}".format(v)

                try:
                    un = sjson["cond_units"][k]
                except KeyError:  # units not defined, or no units for this condition
                    un = ""
                info.append("{0}{1}{2}".format(k, vs, un))
                # Note: should test for filename validity here.
                # See https://stackoverflow.com/questions/9532499/check-whether-a-path-is-valid-in-python-without-creating-a-file-at-the-paths-ta
                # but it looks long. Best is probably to just test write a file
            else:
                if verbose:
                    print(("Warning. {0} not a valid condition".format(k)))
        info = "_".join([_f for _f in info if _f])
    else:
        info = ""

    # ... clean from forbidden characters
    for c in [r"/"]:
        info = info.replace(c, "")

    # ... get extension
    rad, ext = splitext(name)
    if ext == "":
        ext = ".spec"  # default extension

    # ... Write full name
    name = "_".join([_f for _f in [date, rad, info] if _f]) + ext
    fout = join(fold, name)

    # ... Test for existence, replace if needed
    if exists(fout):
        if if_exists_then == "increment":
            if verbose:
                print("Warning. File already exists. Filename is incremented")
            i = 0
            while exists(fout):
                i += 1
                name = "_".join([_f for _f in [date, rad, info, str(i)] if _f]) + ext
                fout = join(fold, name)
        elif if_exists_then == "replace":
            if verbose:
                print(("File exists and will be replaced: {0}".format(name)))
        elif if_exists_then == "ignore":
            if verbose:
                print(("File already exists : {0}. Ignoring".format(name)))
            # main function will return after this function returns.
        else:
            raise ValueError(
                "File already exists {0}. Choose another filename".format(fout)
                + ", or set the `if_exists_then` option to `replace` or ìncrement`"
            )

    return fout


def _compress(s: Spectrum, sjson):
    """removes all quantities that can be regenerated with s.update(), e.g,
    transmittance if abscoeff and path length are given, radiance if emisscoeff
    and abscoeff are given in non-optically thin case, etc.

    Default ``False``
    """

    # TODO: at the moment, a Spectrum with 'radiance_noslit' and 'transmittance_noslit'
    # is saved with 'radiance_noslit' only, if stored with compress, but it is then impossible
    # to reconstruct 'radiance_noslit'.

    from radis.spectrum.rescale import get_redundant

    redundant = get_redundant(s)

    discarded = []
    for key in list(sjson["_q"].keys()):
        if key == "wavespace":
            continue
        if redundant[key]:
            del sjson["_q"][key]
            discarded.append(key)

    if len(discarded) > 0:
        print(
            (
                "Redundant quantities removed: {0}. Use s.update() after ".format(
                    discarded
                )
                + "loading to regenerate them"
            )
        )

    return sjson


# %% Load functions


def load_spec(file, binary=True) -> Spectrum:  # , return_binary_status=False):
    """Loads a .spec file into a :class:`~radis.spectrum.spectrum.Spectrum`
    object. Adds ``file`` in the Spectrum
    :attr:`~radis.spectrum.spectrum.Spectrum.file` attribute.

    Parameters
    ----------
    file: str
        .spec file to load
    binary: boolean
        set to ``True`` if the file is encoded as binary. Default ``True``. Will autodetect
        if it fails, but that may take longer.

    Returns
    -------
    Spectrum : a :class:`~radis.spectrum.spectrum.Spectrum` object

    Examples
    --------
    .. minigallery:: radis.load_spec

    See Also
    --------

    :py:class:`~radis.tools.database.SpecDatabase`,
    :py:meth:`~radis.spectrum.spectrum.Spectrum.store`
    """

    def _load(binary):
        if not binary:
            with open(file, "r") as f:
                sload = json_tricks.load(f, preserve_order=False, ignore_comments=False)
        else:
            with open(file, "rb") as f:
                sload = json_tricks.load(f, preserve_order=False, ignore_comments=False)
        return sload

    # first try to open with given binary info
    try:
        sload = _load(binary)
    # if it fails, retry with different format
    except (TypeError, Exception):
        print(
            (
                "Could not open file {0} with binary={1}. ".format(
                    basename(file), binary
                )
                + "Trying with binary={0}".format(not binary)
            )
        )
        binary = not binary
        sload = _load(binary)
        # if it works:
        print(
            "Worked! Use binary={0} directly in load_spec for faster loading".format(
                binary
            )
        )

    # Test format / correct deprecated format:
    sload, fixed = _fix_format(file, sload)

    # Generate Spectrum
    s = _json_to_spec(sload, file)

    # Auto-update RADIS .spec format
    # ... experimental feature...
    if fixed:
        _update_to_latest_format(s, file, binary)

    return s


def _json_to_spec(sload, file="") -> Spectrum:
    """Builds a Spectrum object from a JSON dictionary. Called by
    :func:`~radis.tools.database.load_spec`.

    Json has been fixed from deprecating changes by _fix_format

    Parameters
    ----------
    sload: dict
        Spectrum object content stored under a dictonary

    Returns
    -------
    Spectrum: a :class:`~radis.spectrum.spectrum.Spectrum` object
    """

    conditions = sload["conditions"]
    # Get quantities
    if "quantities" in sload:
        # old format -saved with tuples (w,I) under 'quantities'): heavier, but
        # easier to generate a spectrum
        quantities = {
            k: (np.array(v[0]), array(v[1])) for (k, v) in sload["quantities"].items()
        }
        warn(
            "File {0}".format(basename(file))
            + " has a deprecrated structure ("
            + "quantities are stored with shared wavespace: uses less space). "
            + "Regenerate database ASAP.",
            DeprecationWarning,
        )
    else:
        quantities = {
            k: (sload["_q"]["wavespace"], v)
            for k, v in sload["_q"].items()
            if k != "wavespace"
        }

    # Generate spectrum:
    waveunit = sload["conditions"]["waveunit"]
    # Only `quantities` and `conditions` is required. The rest is just extra
    # details
    kwargs = {}

    # ... load slit if exists
    if "_slit" in sload:
        slit = {k: np.array(v) for k, v in sload["_slit"].items()}
    else:
        slit = {}

    # ... load lines if exist
    if "lines" in sload:
        df = sload["lines"]
        kwargs["lines"] = df
    else:
        kwargs["lines"] = None

    # ... load populations if exist
    if "populations" in sload:

        # Fix some problems in json-tricks
        # ... cast isotopes to int  (getting sload['populations'] directly doesnt do that)
        kwargs["populations"] = {}
        for molecule, isotopes in sload["populations"].items():
            kwargs["populations"][molecule] = {}
            for isotope, states in isotopes.items():
                try:
                    isotope = int(isotope)
                except ValueError:
                    pass  # keep isotope as it was
                kwargs["populations"][molecule][isotope] = states

    else:
        kwargs["populations"] = None

    # ... load references if exist
    if "references" in sload:
        references = sload["references"]
    else:
        references = {}

    # ... load other properties if exist
    for attr in ["units", "cond_units", "name"]:
        try:
            kwargs[attr] = sload[attr]
        except KeyError:
            kwargs[attr] = None

    s = Spectrum(
        quantities=quantities,
        conditions=conditions,
        wunit=waveunit,
        references=references,
        **kwargs,
    )

    # ... add file
    s.file = basename(file)

    # ... add slit
    s._slit = slit

    return s


def _fix_format(file, sload):
    """Test format / correct deprecated format: The goal is to still be able to
    load old format precomputed spectra, and fix their attribute names.

    Save them again to fix the warnigns definitly.

    Returns
    -------
    json, fixed:  fixed is True if a change was made.
    """

    fixed = False

    # Fix syntax of radis <= 0.1
    # --------------------------

    # ... note: at some points all these tests should simply be removed

    if "q" in sload:
        printr(
            "File {0}".format(basename(file))
            + " has a deprecrated structure (key "
            + "q replaced with _q). Fixed this time, but regenerate "
            + "database ASAP."
        )  # , DeprecationWarning)
        sload["_q"] = sload.pop("q")
        sload["_q_conv"] = sload.pop("q_conv")
        fixed = True

    # Fix _q_conv removal (0.9.30)
    if "_q_conv" in sload:
        if len(sload["_q_conv"]) == 0:
            del sload["_q_conv"]
        else:
            printr(
                "File {0}".format(basename(file))
                + " has a deprecrated structure (key "
                + "_q_conv removed in 0.9.30). Fixed this time, but regenerate "
                + "database for faster loading."
            )
            if not "_q" in sload or len(sload["_q"]) == 0:
                # only convolved quantities; just replace the dict
                sload["_q"] = sload.pop("_q_conv")
            else:
                # here wavespaces arent necessarily the same
                if sload["_q"]["wavespace"] == sload["_q_conv"]["wavespace"]:
                    sload["_q"].update(sload.pop("_q_conv"))
                else:
                    # we need to interpolate
                    raise NotImplementedError

        fixed = True

    try:
        sload["conditions"]["waveunit"]
    except KeyError as err:
        # deprecation: waveunit was named wavespace
        if "wavespace" in sload["conditions"]:
            sload["conditions"]["waveunit"] = sload["conditions"]["wavespace"]
            del sload["conditions"]["wavespace"]
            fixed = True

        else:
            raise KeyError(
                "Spectrum 'conditions' dict should at least have a "
                + "'waveunit' key. Got: {0}".format(list(sload["conditions"].keys()))
            ) from err

    # propagation medium removed in 0.9.22, replaced with 'nm' and 'nm_vac' in
    # waveunit directly
    if "medium" in sload["conditions"]:
        printr(
            "File {0}".format(basename(file))
            + " has a deprecrated structure (key "
            + "medium removed in 0.9.22). Fixing this time, but regenerate "
            + "database ASAP."
        )  # , DeprecationWarning)
        # Fix: rewrite waveunit
        assert "waveunit" in sload["conditions"]
        if sload["conditions"]["waveunit"] == "cm-1":
            pass  # does not change anything, no need to report
        else:  # wavelength is in air or vacuum.
            assert sload["conditions"]["waveunit"] == "nm"
            if sload["conditions"]["medium"] == "air":
                sload["conditions"]["waveunit"] = "nm"
            elif sload["conditions"]["medium"] == "vacuum":
                sload["conditions"]["waveunit"] = "nm_vac"
            else:
                raise ValueError(sload["conditions"]["medium"])
        # fix: delete medium key
        del sload["conditions"]["medium"]
        fixed = True

    if "isotope_identifier" in sload["conditions"]:
        printr(
            "File {0}".format(basename(file))
            + " has a deprecrated structure (key "
            + "isotope_identifier replaced with isotope). Fixed this time, but regenerate "
            + "database ASAP."
        )  # , DeprecationWarning)
        sload["conditions"]["isotope"] = sload["conditions"].pop("isotope_identifier")
        fixed = True

    if "air_pressure_mbar" in sload["conditions"]:
        printr(
            "File {0}".format(basename(file))
            + " has a deprecrated structure (key "
            + "air_pressure_mbar replaced with pressure_mbar). Fixed this time, but regenerate "
            + "database ASAP."
        )  # , DeprecationWarning)
        sload["conditions"]["pressure_mbar"] = sload["conditions"].pop(
            "air_pressure_mbar"
        )
        fixed = True

    if "isotope" in sload["conditions"]:
        isotope = sload["conditions"]["isotope"]
        if not isinstance(isotope, str):
            printr(
                "File {0}".format(basename(file))
                + " has a deprecrated structure (key "
                + "isotope is now a string). Fixed this time, but regenerate "
                + "database ASAP."
            )  # , DeprecationWarning)
            # Fix it:
            sload["conditions"]["isotope"] = ",".join(
                [str(k) for k in list_if_float(isotope)]
            )
            fixed = True

    if "dbpath" in sload["conditions"]:
        dbpath = sload["conditions"]["dbpath"]
        if not isinstance(dbpath, str):
            printr(
                "File {0}".format(basename(file))
                + " has a deprecrated structure (key "
                + "dbpath is now a string). Fixed this time, but regenerate "
                + "database ASAP."
            )  # , DeprecationWarning)
            # Fix it:
            sload["conditions"]["dbpath"] = ",".join(
                [str(k).replace("\\", "/") for k in list_if_float(dbpath)]
            )  # list_if_float or just list??
            fixed = True

    if "selfabsorption" in sload["conditions"]:
        self_absorption = sload["conditions"]["selfabsorption"]
        sload["conditions"]["self_absorption"] = self_absorption
        del sload["conditions"]["selfabsorption"]
        fixed = True

    if "broadening_max_width" in sload["conditions"]:
        broadening_max_width = sload["conditions"]["broadening_max_width"]
        sload["conditions"]["truncation"] = broadening_max_width / 2
        sload["conditions"]["neighbour_lines"] = broadening_max_width / 2
        del sload["conditions"]["broadening_max_width"]
        fixed = True

    # Fix all path names (if / are stored it screws up the JSON loading)
    # -----------------
    def fix_path(key):
        fixed = False
        if key in sload["conditions"]:
            path = sload["conditions"][key]
            if path is not None and not isinstance(path, str):
                printr(
                    "File {0}".format(basename(file))
                    + " has a deprecrated structure (key "
                    + "{0} is now a string). Fixed this time, but regenerate ".format(
                        key
                    )
                    + "database ASAP."
                )  # , DeprecationWarning)
                # Fix it:
                sload["conditions"][key] = path.replace("\\", "/")
                fixed = True
        return fixed

    for param in [
        "database",
        "levelspath",
        "parfuncpath",  # RADIS quantities
        "results_directory",
        "jobName",  # other quantities
    ]:
        fixed += fix_path(param)

    # Fix syntax of radis <= 0.2
    # --------------------------

    # Add 'thermal_equilibrium' in conditions
    if "thermal_equilibrium" not in sload["conditions"]:
        # Hints that the Spectrum has been calculated by a Spectral code:
        if all_in(
            ["Tgas", "Trot", "Tvib", "wstep"], sload["conditions"]
        ) or all_in(  # NEQ / RADIS
            ["Tgas", "Trot", "Tvib", "narray", "jobName"], sload["conditions"]
        ):  # SPECAIR

            def _is_at_equilibrium(conditions):
                # Guess if we're at equilibrium:
                # ... if at Equilibrium, more quantities can be calculated with Kirchoff's law
                # ... using Spectrum.update(). Also, the Spectrum can be stored by deleting
                # ... more redundant parts.
                # ... if not possible to guess, assume False
                try:
                    assert conditions["Tvib"] == conditions["Tgas"]
                    assert conditions["Trot"] == conditions["Tgas"]
                    if "overpopulation" in conditions:
                        assert conditions["overpopulation"] is None
                    assert conditions["self_absorption"]  # is True

                    return True
                except AssertionError:
                    return False
                except KeyError:
                    # Not all keys necessary to decide. Assume False
                    warn(
                        "Missing keys to tell if Spectrum {0} is at equilibrium".format(
                            file
                        )
                        + ". Update spectrum manually"
                    )
                    return None

            equilibrium = _is_at_equilibrium(sload["conditions"])

            if equilibrium is not None:
                sload["conditions"]["thermal_equilibrium"] = equilibrium
                fixed = True
                printr(
                    "File {0}".format(basename(file))
                    + " has a deprecrated structure ("
                    + "thermal_equilibrium not defined). Fixed it this time (guessed {0})".format(
                        equilibrium
                    )
                    + ", but regenerate file ASAP."
                )  # , DeprecationWarning)

    # Fix lines format HITRAN_CLASS_1 molecules
    if "lines" in sload and sload["lines"] is not None:
        lines = sload["lines"]
        from radis.db.classes import HITRAN_CLASS1, get_molecule

        if "v1u" in lines and get_molecule(lines.id.iloc[0]) in HITRAN_CLASS1:
            printr(
                "File {0}".format(basename(file))
                + " has a deprecrated structure "
                + "(v1u in lines is now called vu). Fixed this time, but regenerate "
                + "database ASAP."
            )
            # Fix it:
            lines.rename(columns={"v1u": "vu", "v1l": "vl"}, inplace=True)

    # Fix syntax of RAdis <= 0.9.26
    # -----------------------------

    # Fix adimensioned units
    if "units" in sload:
        for var, unit in sload["units"].items():
            if unit in ["I/I0", "-ln(I/I0)", "eps"]:
                printr(
                    "File {0}".format(basename(file))
                    + " has a deprecrated structure "
                    + "(adimensioned units are now stored as ''). Fixed this time, but regenerate "
                    + "database ASAP."
                )
                sload["units"][var] = ""
            if "cm_1" in unit:
                printr(
                    "File {0}".format(basename(file))
                    + " has a deprecrated structure "
                    + "(cm_1 is now written cm-1''). Fixed this time, but regenerate "
                    + "database ASAP."
                )
                sload["units"][var] = sload["units"][var].replace("cm_1", "cm-1")

    return sload, fixed


def _update_to_latest_format(s, file, binary):
    """experimental feature Used to autoupdate .spec files to the latest
    format, by simply saving them again once they're loaded and fixed. Warning!
    Better have a copy of your files before that, or a way to regenerate them.

    Parameters
    ----------
    s: a Spectrum object
        already loaded, already updated. See fixes in :func:`~radis.tools.database._fix_format`
    file: str
        current storage file of the object
    binary: boolean
        the binary format that was used to open the object. Keep the same.
        if ``True``, saves under binary format. Faster and takes less space.
        If ``2``.

    Examples
    --------
    add to the start of your script::

        import radis
        radis.config["AUTO_UPDATE_SPEC"] = True
    """

    import radis

    if radis.config["AUTO_UPDATE_SPEC"]:

        assert exists(file)

        s.store(file, compress=binary, if_exists_then="replace", discard=[])

        printr("File {0} auto-updated to latest RADIS format".format(file))

        return


def plot_spec(file, what="radiance", title=True, **kwargs):
    """Plot a .spec file. Uses the
    :py:meth:`~radis.spectrum.spectrum.Spectrum.plot` method internally.

    Parameters
    ----------
    file: str,  or Spectrum object
        .spec file to load, or Spectrum object directly

    Other Parameters
    ----------------
    kwargs: dict
        arguments forwarded to :meth:`~radis.spectrum.spectrum.Spectrum.plot`


    Returns
    -------
    fig: matplotlib figure
        where the Spectrum has been plotted

    See Also
    --------

    :py:meth:`~radis.spectrum.spectrum.Spectrum.plot`
    """

    import matplotlib.pyplot as plt

    if isinstance(file, str):
        s = load_spec(file)
    elif isinstance(file, Spectrum):
        s = file
    else:
        raise ValueError(
            "file should be a string, or a Spectrum object directly. "
            + "Got {0}".format(type(file))
        )

    try:
        s.plot(what, **kwargs)
    except KeyError:
        try:
            print((sys.exc_info()[0], sys.exc_info()[1]))
            s.plot(what + "_noslit", **kwargs)  # who knows maybe it will work :)
            print(("Printing {0} instead".format(what + "_noslit")))
        except:
            print((sys.exc_info()[0], sys.exc_info()[1]))
            # Plot something
            s.plot(s.get_vars()[0], **kwargs)
            print(("Printing {0} instead".format(s.get_vars()[0])))

    if title and s.file:
        plt.title(basename(s.file))
        plt.tight_layout()

    return plt.gcf()


# %% SpecList class
# ... loads a list of spectra and manipulate them


class SpecList(object):
    def __init__(self, *spectra, **kwargs):
        """A list of Spectrum, with various methods to manage them.

        .. warning::     still new in 0.9.17
        """

        # get defaults
        verbose = kwargs.pop("verbose", True)

        self.df = self._create_df(*spectra)
        self.verbose = verbose

    def _create_df(self, *spectra):
        """Returns Dataframe of spectra with conditions.

        Parameters
        ----------

        spectra: list of Spectrum
            all spectra
        """

        # TODO: deal with case where SpecList(list) was called (see below for first attempt)

        #        # in case SpecList(list) was called instead of SpecList(*list)
        #        if isinstance(spectra, tuple):
        #            assert len(spectra) == 1
        #            spectra = spectra[0]

        db = []

        for s in spectra:
            params = s.get_conditions().copy()
            params.update({"id": id(s), "Spectrum": s})
            db.append(params)

        return pd.DataFrame(db)

    def conditions(self):
        """Show conditions in database."""

        cond = list(self.df.columns)

        if len(cond) > 0:
            if "file" in cond:
                cond.remove("file")
            cond.remove("Spectrum")

        return cond

    def see(self, columns=None, *args):
        """Shows Spectrum database with all conditions (``columns=None``) or
        specific conditions.

        Parameters
        ----------

        columns: str, list of str, or None
            shows the conditions value for all cases in database. If None, all
            conditions are shown. Default ``None``
            e.g.::

                db.see(['Tvib', 'Trot'])

        Notes
        -----

        Makes the 'file' column the index, and also discard the 'Spectrum' column
        (that holds all the data) for readibility
        """

        if len(self) == 0:
            raise ValueError("Database is empty")

        if isinstance(columns, str):
            columns = [columns] + [k for k in args]

        dg = self.df.set_index("file")  # note that this is a copy already.
        # dont try to modify the output of "see"
        del dg["Spectrum"]  # for visibility"

        if columns is None:
            return dg

        for c in columns:
            if not c in self.df.columns:
                raise ValueError(
                    "`{0}` is not a column name. Use one of {1}".format(
                        c, self.df.columns
                    )
                )

        return dg.reindex(columns=columns)

    def view(self, columns=None, *args):
        """alias of :py:meth:`~radis.tools.database.SpecList.see`

        See Also
        --------

        :py:meth:`~radis.tools.database.SpecList.see`
        """

        return self.see(columns=columns, *args)

    def map(self, function):
        """Apply ``function`` to all Spectra in database.

        Examples
        --------

        Add a missing parameter::

            db = SpecDatabase('...')

            def add_condition(s):
                s.conditions['exp_run'] = 1
                return s

            db.map(add_condition)

        .. note::

            spectra are not changed on disk. If you want to update on disk
            you may want to combine map() followed by :py:meth:`~radis.tools.database.SpecDatabase.compress_to`

        Example ::

            # See length of all spectra :

            db.map(lambda s: print(len(s)))


            # Resample all on spectrum of minimum wstep
            s_wstep_min = db.get(wstep=float(db.see("wstep").min()))[0]

            db.map(lambda s: s.resample(s_wstep_min))

            # Export to a new database:
            db.compress_to(db.path+'_interp')

        """
        for s in self:
            function(s)

    def _load_existing_files(self, files):
        """Loadspectra already registered in the database.

        Parameters
        ----------
        files : list

        Returns
        -------
        list
            of Spectrum

        """
        # Parallel loading
        if len(files) > self.minimum_nfiles:
            nJobs = self.nJobs
            batch_size = self.batch_size
            if self.verbose:
                print(
                    "*** Loading the database with {0} ".format(
                        nJobs if nJobs > 0 else f"all minus {abs(nJobs+1)}"
                    )
                    + "processor(s) ({0} files)***".format(len(files))
                )

            def funLoad(f):
                spectrum = load_spec(join(self.path, f), binary=self.binary)
                indx = self.df.index[self.df.index[self.df["file"] == f]].to_list()
                if len(indx) > 1:
                    raise ValueError(
                        "Multiple files found for the name {} of indexes {}. Update the database manually".format(
                            f, indx
                        )
                    )
                # self.df.loc[indx[0],"Spectrum"] = spectrum
                return indx[0], spectrum

            spec_loaded = Parallel(
                n_jobs=nJobs, batch_size=batch_size, verbose=5 * self.verbose
            )(delayed(funLoad)(f) for f in files)
            for val in spec_loaded:
                self.df.loc[val[0], "Spectrum"] = val[1]
        else:
            for f in files:
                idx = self.df.index[self.df.index[self.df["file"] == f]].to_list()
                if len(idx) > 1:
                    raise ValueError(
                        "Multiple files found for the name {} of indexes {}. Update the database manually".format(
                            f, idx
                        )
                    )
                spectrum = load_spec(join(self.path, f), binary=self.binary)
                self.df.loc[idx[0], "Spectrum"] = spectrum
        return list(self.df["Spectrum"])

    def get(self, conditions="", **kwconditions):
        """Returns a list of spectra that match given conditions.

        Parameters
        ----------
        database: list of Spectrum objects
            the database
        conditions: str
            a list of conditions. Example::

                db.get('Tvib==3000 & Trot==1500')

        kwconditions: dict
            an unfolded dict of conditions. Example::

                db.get(Tvib=3000, Trot=1500)

        Other Parameters
        ----------------
        inplace: ``bool``
            if True, return the actual object in the database. Else, return
            copies. Default ``False``
        verbose: ``bool``
            more blabla
        scale_if_possible: ``bool``
            if ``True``, spectrum is scaled for parameters that can be computed
            directly from spectroscopic quantities (e.g: ``'path_length'``,
            ``'molar_fraction'``). Default ``False``

        Returns
        -------
        out: list of Spectrum

        Examples
        --------
        ::

            spec_list = db.get('Tvib==3000 & Trot==1300')

        or::

            spec_list = db.get(Tvib=3000, Trot=1300)


        See Also
        --------
        :meth:`~radis.tools.database.SpecList.get_unique`,
        :meth:`~radis.tools.database.SpecList.get_closest`,
        ;py:meth:`~radis.tools.database.SpecDatabase.interpolate`,
        :meth:`~radis.tools.database.SpecList.items`
        """

        ## type: bool, default False
        kwconditions.pop("verbose", True)  #: type: bool
        inplace = kwconditions.pop("inplace", False)
        scale_if_possible = kwconditions.pop("scale_if_possible", False)

        # Test inputs
        for (k, _) in kwconditions.items():
            if not k in self.df.columns:
                raise ValueError(
                    "{0} not a correct condition name. Use one of: {1}".format(
                        k, self.df.columns
                    )
                )
        if len(self.df) == 0:
            warn("Empty database")
            return []

        if scale_if_possible:
            scaled_inputs = {}
            # Remove scalable inputs from required keys (we scale them at the end instead)
            for k in _scalable_inputs:
                if k in kwconditions and len(kwconditions) > 1:
                    scaled_inputs[k] = kwconditions.pop(k)

        # Unique condition method
        if conditions != "" and kwconditions != {}:
            raise ValueError(
                "Please choose one of the two input format (str or dict) exclusively"
            )

        if conditions == "" and kwconditions == {}:  # get all spectra
            # Get all unloaded Spectrum objects and load them
            files = self.df["file"][self.df["Spectrum"].isnull()]
            self._load_existing_files(files)
            out = list(self.df["Spectrum"])

        else:
            # Find Spectrum that match conditions
            if conditions != "":  # ... with input conditions query directly
                dg = self.df.query(conditions)
            else:  # ... first write input conditions query
                query = []
                for (k, v) in kwconditions.items():
                    if isinstance(v, str):
                        query.append("{0} == r'{1}'".format(k, v))
                    elif v is None:
                        # query "k == None" doesn't work. We use a workaround,
                        # checking if the column is different from itself (i.e. : is None):
                        # https://stackoverflow.com/a/32207819/5622825
                        query.append(f"{k} != {k}")
                    else:
                        #                    query.append('{0} == {1}'.format(k,v))
                        query.append("{0} == {1}".format(k, v.__repr__()))
                        # ... for som reason {1}.format() would remove some digit
                        # ... to floats in Python2. Calling .__repr__() keeps
                        # ... the correct format, and has no other consequences as far
                        # ... as I can tell

                # There is a limitation in numpy: a max of 32 arguments is required.
                # Below we write a workaround when the Spectrum has more than 32 conditions
                if len(query) < 32:
                    query = " & ".join(query)
                    if __debug__:
                        printdbg("Database query: {0}".format(query))
                    dg = self.df.query(query)
                else:
                    # cut in <32-long parts
                    N = len(query) // 32 + 1
                    querypart = " & ".join(query[::N])
                    dg = self.df.query(querypart)
                    for i in range(1, N + 1):
                        querypart = " & ".join(query[i::N])
                        if __debug__:
                            printdbg("Database query: {0}".format(querypart))
                        dg = dg.query(querypart)
            # Get all unloaded Spectrum objects and load them
            files = dg["file"][dg["Spectrum"].isnull()]

            self._load_existing_files(files)
            out = list(self.df.loc[dg.index, "Spectrum"])

        if not inplace:
            out = [s.copy() for s in out]

        # Scale scalable conditions
        if scale_if_possible:
            for s in out:
                for k, v in scaled_inputs.items():
                    if k == "path_length":
                        s.rescale_path_length(v)
                    elif k == "mole_fraction":
                        s.rescale_mole_fraction(v)
                    else:
                        raise KeyError("cant rescale this: {0}".format(k))

        return out

    def get_unique(self, conditions="", scale_if_possible=False, **kwconditions):
        """Returns a spectrum that match given conditions.

        Raises an error if the spectrum is not unique.

        Parameters
        ----------
        args:
            see :py:meth:`~radis.tools.database.SpecList.get` for more details

        Returns
        -------
        s: Spectrum

        See Also
        --------
        :py:meth:`~radis.tools.database.SpecList.get`,
        :py:meth:`~radis.tools.database.SpecList.get_closest`,
        ;py:meth:`~radis.tools.database.SpecDatabase.interpolate`
        """

        out = self.get(conditions, scale_if_possible=scale_if_possible, **kwconditions)

        if len(out) == 0:
            # Give a better error message before crashing:
            try:
                prevVerbose = kwconditions.get("verbose", None)
                kwconditions["verbose"] = True
                self.get_closest(**kwconditions)  # note: wont work with conditions=..
            except:
                pass
            finally:
                if prevVerbose is not None:
                    kwconditions["verbose"] = prevVerbose
                raise ValueError(
                    "Spectrum not found. See closest above. Use db.get_closest(). You could also try db.interpolate()"
                )
        elif len(out) > 1:
            raise ValueError(
                "Spectrum is not unique ({0} match found)".format(len(out))
                + ' Think about using "db.find_duplicates()"'
            )
        else:
            return out[0]

    def get_closest(self, scale_if_possible=True, **kwconditions):
        """Returns the Spectra in the database that is the closest to the input
        conditions.

        Note that for non-numeric values only equals should be given.
        To calculate the distance all numeric values are scaled by their
        mean value in the database

        Parameters
        ----------
        kwconditions: named arguments
            i.e: ``Tgas=300, path_length=1.5``
        scale_if_possible: boolean
            if ``True``, spectrum is scaled for parameters that can be computed
            directly from spectroscopic quantities (e.g: ``'path_length'``, ``'molar_fraction'``).
            Default ``True``

        Other Parameters
        ----------------
        verbose: boolean
            print messages. Default ``True``
        inplace: boolean
            if ``True``, returns the actual object in database. Else, return a copy.
            Default ``False``

        See Also
        --------
        :meth:`~radis.tools.database.SpecList.get`,
        :meth:`~radis.tools.database.SpecList.get_unique`,
        ;py:meth:`~radis.tools.database.SpecDatabase.interpolate`
        """
        #        split_columns: list of str.
        #            slits a comma separated column in multiple columns, and number them.
        #            Typically::
        #
        #                db.get_closest(..., split_columns=['Tvib'])    # in multi temperature modes
        #
        #            splits ``Tvib:'1200,1300,1000'``, in ``Tvib1:1200, Tvib2:1300, Tvib3:1000``
        #            Default ``[]``

        # TODO: make it possible to choose only certain parameters, fix others.
        # Maybe first generate a SpecDatabase of output of spectra with get(),
        # then get_closest() for the parameters we dont want to fix

        # Test inputs
        if kwconditions == {}:
            raise ValueError("Please specify filtering conditions. e.g: Tgas=300")

        # Inputs:
        verbose = kwconditions.pop("verbose", True)  #: type: bool
        inplace = kwconditions.pop("inplace", False)  #: type: bool
        #        split_columns = kwconditions.pop('split_columns', [])     #: type: bool

        # Check all conditions exist
        for (k, _) in kwconditions.items():
            if not k in self.df.columns:
                raise ValueError(
                    "{0} not a correct condition name. Use one of: {1}".format(
                        k, self.df.columns
                    )
                )

        dg = self.df.reindex(columns=[k for k in self.df.columns if k != "Spectrum"])

        if scale_if_possible:
            # Remove scalable inputs from distance calculation variables (unless
            # they're the last variable, because then it may screw things up)
            for k in _scalable_inputs:
                if len(dg.columns) == 1:
                    break
                try:
                    del dg[k]
                except KeyError:
                    pass

        mean = dict(dg.mean())
        #        std = dict(dg.std())
        assert "_d" not in dg.columns
        dg["_d"] = 0  # distance (squared, actually)
        try:
            for k, v in kwconditions.items():
                if not k in dg.columns:
                    continue
                #            # add distance to all conditions planes. We regularize the different
                #            # dimensions by working on normalized quantities:
                #            #     a  ->   (a-mean)/std  € [0-1]
                #            # Distance becomes:
                #            #     d^2 ->  sum((a-target)/std)^2
                #            # Problem when std == 0! That means this dimension is not discrimant
                #            # anyway
                #            if std[k] == 0:
                #                # for this conditions all parameters have the same value.
                #                dg['_d'] += (dg[k]-v)**2
                #            else:
                #                dg['_d'] += (dg[k]-v)**2/(std[k])**2
                # Eventually I chose to scale database with mean only (there was
                # an obvious problem with standard deviation scaling in the case of
                # a non important feature containing very close datapoints that would
                # result in inappropriately high weights)

                try:
                    dg["_d"] += (dg[k] - v) ** 2 / mean[k] ** 2
                except TypeError as err:
                    # Deal with case where formats dont match:
                    try:
                        dg[k] - v
                    except TypeError as err2:
                        print(sys.exc_info())
                        raise TypeError(
                            "An error occured (see above) when calculating "
                            + f"(dg[{k}] - {v}). Example: "
                            + f"({dg[k].iloc[0]} - {v}). "
                            + "Check that your requested conditions match "
                            + "the database format"
                        ) from err2
                    else:
                        raise err

            # Get spectrum with minimum distance to target conditions
            ## type: Spectrum
            sout = self.df.loc[dg["_d"].idxmin(), "Spectrum"]
            # Load the spectrum if not already done
            if sout is None:
                spectrum = load_spec(
                    join(self.path, self.df.loc[dg["_d"].idxmin(), "file"]),
                    binary=self.binary,
                )
                self.df.loc[dg["_d"].idxmin(), "Spectrum"] = spectrum
                sout = spectrum
            # Note @EP 07/12/20 : do we have the same index as dg ?? (created with reindex on L1335)  #  TODO
        finally:
            del dg["_d"]

        if not inplace:
            sout = sout.copy()

        # Scale scalable conditions
        if scale_if_possible:
            if "path_length" in kwconditions:
                sout.rescale_path_length(kwconditions["path_length"])
            if "mole_fraction" in kwconditions:
                sout.rescale_mole_fraction(kwconditions["mole_fraction"])

        if verbose > 1:
            print(
                (
                    "{:^7} \t".format(self.df.molecule[0])
                    + "\t".join(["{0}".format(k) for k in kwconditions.keys()])
                )
            )
            print(
                (
                    "Look up \t"
                    + "\t".join(["{0:.3g}".format(v) for v in kwconditions.values()])
                )
            )
            print(
                (
                    "Got     \t"
                    + "\t".join(
                        [
                            "{0:.3g}".format(sout.conditions[k])
                            for k in kwconditions.keys()
                        ]
                    )
                )
            )

        return sout

    def get_items(self, condition):
        """Returns all Spectra in database under a dictionary; indexed by
        ``condition``

        Requires that ``condition`` is unique

        Parameters
        ----------
        condition: str
            condition. Ex: ``Trot``

        Returns
        -------
        out: dict
            {condition:Spectrum}

        Examples
        --------
        ::

            db.get_items("Tgas")

        See Also
        --------
        :meth:`~radis.tools.database.SpecList.to_dict`,
        :meth:`~radis.tools.database.SpecList.get`,
        :py:meth:`~radis.tools.database.SpecList.create_fname_grid`
        """

        if not self.df[condition].is_unique:
            raise ValueError(
                "Values in {0} must be unique to use get_items(). Got {1}".format(
                    condition, self.see(condition)[self.see(condition).duplicated()]
                )
            )

        return dict(zip(self.df[condition], self.df.Spectrum))

    def create_fname_grid(self, conditions):
        """Create a 2D-grid of filenames for the list of parameters ``conditions``

        Examples
        --------
        ::

            db.create_fname_grid(["Tgas", "pressure_mbar"])

        See Also
        --------
        :py:meth:`~radis.tools.database.SpecList.get_items`

        """

        gridsize = [len(self.df[cond].unique()) for cond in conditions]

        return (
            self.see(conditions)
            .reset_index()
            .set_index(conditions)
            .values.reshape(gridsize)
        )

    def __iter__(self):
        """Iterate over all Spectra in database.

        .. warning::

            returns the inplace object directly. If you modify them, the Spectra
            are modified

        Examples
        --------
        Print name of all Spectra in dictionary::

            db = SpecDatabase('.')
            for s in db:
                print(s.name)

        Note that this is equivalent to::

            db = SpecDatabase('.')
            for s in db.get():
                print(s.name)

        or::

            db = SpecDatabase('.')
            db.map(lambda s: print(s.name))

        See Also
        --------
        :meth:`~radis.tools.database.SpecList.keys`,
        :meth:`~radis.tools.database.SpecList.values`,
        :meth:`~radis.tools.database.SpecList.items`,
        :meth:`~radis.tools.database.SpecList.to_dict`
        """

        return self.get(inplace=True).__iter__()

    def keys(self):
        """Iterate over all {path} in database.

        See Also
        --------

        :meth:`~radis.tools.database.SpecList.values`,
        :meth:`~radis.tools.database.SpecList.items`,
        :meth:`~radis.tools.database.SpecList.to_dict`
        """

        return list(self.to_dict().keys())

    def values(self):
        """Iterate over all {Spectrum} in database.

        See Also
        --------

        :meth:`~radis.tools.database.SpecList.keys`,
        :meth:`~radis.tools.database.SpecList.items`,
        :meth:`~radis.tools.database.SpecList.to_dict`
        """

        return list(self.to_dict().values())

    def items(self):
        """Iterate over all :py:class:`~radis.spectrum.spectrum.Spectrum` in
        database.

        Examples
        --------
        Print name of all Spectra in dictionary::

            db = SpecDatabase('.')
            for path, s in db.items():
                print(path, s.name)

        Update all spectra in current folder with a new condition ('author')::

            db = SpecDatabase('.')
            for path, s in db.items():
                s.conditions['author'] = 'me'
                s.store(path, if_exists_then='replace')

        See Also
        --------
        :meth:`~radis.tools.database.SpecList.keys`,
        :meth:`~radis.tools.database.SpecList.values`,
        :meth:`~radis.tools.database.SpecList.to_dict`
        """

        return list(self.to_dict().items())

    def plot(self, nfig=None, legend=True, **kwargs):
        """Plot all spectra in database.

        Parameters
        ----------
        nfig: str, or int, or ``None``
            figure to plot on. Default ``None``: creates one

        Other Parameters
        ----------------
        kwargs: dict
            parameters forwarded to the Spectrum :meth:`~radis.spectrum.spectrum.plot`
            method
        legend: bool
            if ``True``, plot legend.

        Returns
        -------
        fig, ax: matplotlib figure and ax
            figure

        Examples
        --------
        Plot all spectra in a folder::

            db = SpecDatabase('my_folder')
            db.plot(wunit='nm')

        See Also
        --------
        Spectrum :py:meth:`~radis.spectrum.spectrum.Spectrum.plot` method
        """

        import matplotlib.pyplot as plt

        fig = plt.figure(num=nfig)
        ax = fig.gca()
        for s in self:
            s.plot(nfig="same", **kwargs)

        if legend:
            plt.legend(loc="best")

        return fig, ax

    def plot_cond(self, cond_x, cond_y, z_value=None, nfig=None):
        """Plot database conditions available:

        Parameters
        ----------
        cond_x, cond_y: str
            columns (conditions) of database.
        z_value: array, or None
            if not None, colors the 2D map with z_value. z_value is ordered
            so that z_value[i] corresponds to row[i] in database.

        Examples
        --------
        ::
            >>> db.plot(Tvib, Trot)     # plot all points calculated


            >>> db.plot(Tvib, Trot, residual)     # where residual is calculated by a fitting
                                                  # procedure...

        .. minigallery:: radis.tools.database.SpecList.plot_cond
        """
        # %%

        import matplotlib.pyplot as plt

        from radis.misc.plot import fix_style, set_style

        x = self.df[cond_x]
        y = self.df[cond_y]

        # Default
        set_style()
        fig = plt.figure(num=nfig)
        ax = fig.gca()
        ax.plot(x, y, "ok")
        ax.set_xlabel(cond_x)
        ax.set_ylabel(cond_y)
        title = self.name

        # Overlay color
        if z_value is not None:
            if type(z_value) is str:
                z = self.df[z_value]
            else:
                z = z_value

            assert len(z) == len(self.df)

            z = np.array(z) ** 0.5  # because the lower the better
            #            norm = cm.colors.Normalize(vmax=z.max(), vmin=z.min())
            #            cmap = cm.PRGn

            xarr = np.linspace(min(x), max(x))
            yarr = np.linspace(min(y), max(y))
            mx, my = np.meshgrid(xarr, yarr)
            zgrid = griddata((x, y), z, (mx, my), method="linear", fill_value=np.nan)
            levels = np.linspace(min(z), max(z), 20)
            ax.contourf(
                mx,
                my,
                zgrid,
                levels=levels,
                # linewidths=1,linestyles='dashed',
                extend="both",
            )

            ix0, iy0 = np.where(zgrid == zgrid.min())
            plt.plot(
                (xarr[ix0], xarr[ix0]), (yarr.min(), yarr.max()), color="white", lw=0.5
            )
            plt.plot(
                (xarr.min(), xarr.max()), (yarr[iy0], yarr[iy0]), color="white", lw=0.5
            )
            # lbls = plt.clabel(cs0, inline=1, fontsize=16,fmt='%.0fK',colors='k',manual=False)     # show labels

        plt.title(title)
        plt.tight_layout()
        fix_style()

    def __len__(self):
        return len(self.df)

    def __str__(self):
        return self.see().__str__()


# %% Database class
# ... loads database and manipulate it
# ... similar to SpecList, but associated and synchronized with a folder
# ... on the disk


class SpecDatabase(SpecList):
    """A Spectrum Database class to manage them all.

    It basically manages a list of Spectrum JSON files, adding a Pandas
    dataframe structure on top to serve as an efficient index to visualize
    the spectra input conditions, and slice through the Dataframe with
    easy queries

    Similar to :class:`~radis.tools.database.SpecList`, but associated and
    synchronized with a folder

    Parameters
    ----------
    path: str
        a folder to initialize the database
    filt: str
        only consider files ending with ``filt``. Default ``.spec``
    binary: boolean
        if ``True``, open Spectrum files as binary files. If ``False`` and it fails,
        try as binary file anyway. Default ``False``.
    lazy_loading: bool``
        If ``True``, load only the data from the summary csv file and the spectra will
        be loaded when accessed by the get functions. If ``False``, load all
        the spectrum files. If ``True`` and the summary .csv file does not exist,
        load all spectra

    Other Parameters
    ----------------
    *input for :class:`~joblib.parallel.Parallel` loading of database*

    nJobs: int
        Number of processors to use to load a database (usefull for big
        databases). BE CAREFULL, no check is done on processor use prior
        to the execution ! Default ``-2``: use all but 1 processors.
        Use ``1`` for single processor.
    batch_size: int or ``'auto'``
        The number of atomic tasks to dispatch at once to each
        worker. When individual evaluations are very fast, dispatching
        calls to workers can be slower than sequential computation because
        of the overhead. Batching fast computations together can mitigate
        this. Default: ``'auto'``

    More information in :class:`joblib.parallel.Parallel`

    Examples
    --------

    ::

        >>> db = SpecDatabase(r"path/to/database")     # create or loads database

        >>> db.update()  # in case something changed
        >>> db.see(['Tvib', 'Trot'])   # nice print in console

        >>> s = db.get('Tvib==3000')[0]  # get a Spectrum back
        >>> db.add(s)  # update database (and raise error because duplicate!)

    Note that :py:class:`~radis.lbl.factory.SpectrumFactory` can be configured to
    automatically look-up and update a database when spectra are calculated.

    The function to auto retrieve a Spectrum from database on calculation
    time is a method of DatabankLoader class

    You can see more examples on the :ref:`Spectrum Database section <label_spectrum_database>`
    of the website.

    .. notes:
        Database interaction is based on Pandas query functions. As such, it
        requires all conditions to be either float, string, or boolean. List
        won't work!

    See Also
    --------
    :func:`~radis.tools.database.load_spec`,
    :meth:`~radis.spectrum.spectrum.Spectrum.store`

    Methods to retrieve objects:

    :meth:`~radis.tools.database.SpecList.get`,
    :meth:`~radis.tools.database.SpecList.get_closest`,
    :meth:`~radis.tools.database.SpecList.get_unique`,

    Methods to manipulate the SpecDatabase:

    :meth:`~radis.tools.database.SpecList.see`,
    :meth:`~radis.tools.database.SpecDatabase.update`,
    :meth:`~radis.tools.database.SpecDatabase.add`,
    :meth:`~radis.tools.database.SpecDatabase.compress_to`,
    :meth:`~radis.tools.database.SpecDatabase.find_duplicates`

    Compare another Spectrum to all spectra in the database:

    :meth:`~radis.tools.database.SpecDatabase.fit_spectrum`


    .. minigallery:: radis.SpecDatabase
    """

    def __init__(
        self,
        path=".",
        filt=".spec",
        add_info=None,
        add_date="%Y%m%d",
        verbose=True,
        binary=True,
        nJobs=-2,
        batch_size="auto",
        lazy_loading=True,
        update_register_only=False,
    ):
        # TODO @devs: generate a SpecDatabase from a dict.
        # use the key of the dict insted of the file.

        # Assert name looks like a directory
        name, ext = splitext(str(path))

        if ext != "":
            raise ValueError("Database should be a directory: {0}".format(path))

        self.name = basename(abspath(name))
        self.path = path

        if not exists(path):
            # create it
            os.mkdir(path)
            if verbose:
                print(
                    ("Database {0} initialised in {1}".format(self.name, dirname(path)))
                )
        else:
            if verbose:
                print(("Loading database {0}".format(self.name)))

        #        self.df = None               # created in SpecList.__init__()
        #        self.verbose = verbose       # created in SpecList.__init__()
        self.binary = binary

        # default
        self.add_info = add_info
        self.add_date = add_date

        # loading parameters
        self.nJobs = nJobs
        self.batch_size = batch_size
        self.minimum_nfiles = (
            10  #: type: int. If there are less files, don't use parallel mode.
        )

        # init with 0 spectra
        super(SpecDatabase, self).__init__()
        # now load from the folder with the update() function
        if not lazy_loading:
            return self.update(
                force_reload=True, filt=filt, update_register_only=update_register_only
            )
        else:
            if verbose >= 2:
                print(
                    "Database initialized without loading spectra. They will be loaded when accessed by the get function"
                )
            csv_file = [f for f in os.listdir(self.path) if f.endswith(".csv")]
            if len(csv_file) == 0:
                # load from the folder and initialize the csv with the update() function
                return self.update(force_reload=True, filt=filt)
            if len(csv_file) > 1:
                raise ValueError(
                    "Only 1 csv file should be in the database directory but {0} found : {1}".format(
                        len(csv_file), csv_file
                    )
                    + ". Clean the files then initialize the database with lazy_loading=False"
                )
            spec_files = [f for f in os.listdir(self.path) if f.endswith(".spec")]
            # Initialize database without loading the spectra
            self.df = read_conditions_file(
                join(self.path, csv_file[0]), verbose=self.verbose
            )

            if len(self.df["file"]) != len(spec_files):
                raise ValueError(
                    "Number of files stored in {0} ({1}) doesn't match the number of files found in the database folder ({2}). Update the csv first by setting lazy_loading=False".format(
                        csv_file[0], len(self.df["file"]), len(spec_files)
                    )
                )
            # Check file was not modified since it was registered in the .csv
            for f in spec_files:
                if not f in self.df["file"].values:
                    raise ValueError(
                        "{} not found in the csv registry.  Update the csv first with `SpecDatabase(..., lazy_loading=False)`".format(
                            f
                        )
                    )
                else:
                    if "last_modified" in self.df:
                        if round(
                            self.df["last_modified"][
                                self.df.index[self.df["file"] == f].to_list()[0]
                            ],
                            2,
                        ) != round(
                            os.path.getmtime(join(self.path, f)), 2
                        ):  # take rounded numbers to prevent errors
                            raise ValueError(
                                "Spectrum {} last modification don't match the date stored in the csv. Update the database first with `SpecDatabase(..., lazy_loading=False)`".format(
                                    f
                                )
                            )
                    else:
                        raise ValueError(
                            "Csv file wasn't created with last modifications informations. Update the database first with `SpecDatabase(..., lazy_loading=False)`."
                        )

    def update(self, force_reload=False, filt=".spec", update_register_only=False):
        """Reloads database, updates internal index structure and export it in
        ``<database>.csv``.

        Parameters
        ----------
        force_reload: boolean
            if ``True``, reloads files already in database. Default ``False``
        filt: str
            only consider files ending with ``filt``. Default ``.spec``

        Other Parameters
        ----------------
        update_register_only: bool
            if ``True``, load files and update csv but do not keep the Spectrum in memory.
            Default ``False``

        Notes
        -----
        Can be loaded in parallel using joblib by setting the `nJobs` and `batch_size`
        attributes of :class:`~radis.tools.database.SpecDatabase`. See :class:`joblib.parallel.Parallel`
        for information on the arguments
        """

        path = self.path

        if force_reload:
            # Reloads whole database  (necessary on database init to create self.df
            files = [join(path, f) for f in os.listdir(path) if f.endswith(filt)]
            self.df = self._load_new_files(
                files=files, update_register_only=update_register_only
            )
        else:
            dbfiles = list(self.df["file"])
            files = [
                join(path, f)
                for f in os.listdir(path)
                if f not in dbfiles and f.endswith(filt)
            ]
            # no parallelization here because the number of files is supposed to be small
            for f in files:
                self.df = self.df.append(
                    self._load_new_file(
                        f, binary=self.binary, update_register_only=update_register_only
                    ),
                    ignore_index=True,
                )

        # Print index
        self.print_index()

    def compress_to(self, new_folder, compress=True, if_exists_then="error"):
        """Saves the Database in a new folder with all Spectrum objects under
        compressed (binary) format. Read/write is much faster. After the
        operation, a new database should be initialized in the new_folder to
        access the new Spectrum.

        Parameters
        ----------
        new_folder: str
            folder where to store the compressed SpecDatabase. If doesn't exist,
            it is created.
        compress: boolean, or 2
            if ``True``, saves under binary format. Faster and takes less space.
            If ``2``, additionaly remove all redundant quantities.
        if_exists_then: ``'increment'``, ``'replace'``, ``'error'``, ``'ignore'``
            what to do if file already exists. If ``'increment'`` an incremental digit
            is added. If ``'replace'`` file is replaced (!). If ``'ignore'`` the
            Spectrum is not added to the database and no file is created.
            If ``'error'`` (or anything else) an error is raised. Default ``'error'``.

        See Also
        --------
        :meth:`~radis.tools.database.SpecDatabase.find_duplicates`
        """

        # dont allow to compress to the same folder (doesnt make sense, and will
        # create conflicts with same name files)
        assert abspath(self.path) != abspath(new_folder)

        if not exists(new_folder):
            os.makedirs(new_folder)

        for file, s in self.items():  # loop over all spectra
            s.store(
                join(new_folder, file), compress=compress, if_exists_then=if_exists_then
            )

        if self.verbose:
            print("Database compressed to {0}".format(new_folder))

    def print_index(self, file=None):
        if file is None:
            file = join(self.path, self.name + ".csv")

        if len(self) > 0:
            try:
                self.see().to_csv(file)
            except PermissionError:
                warn(
                    "Database index could not be updated: {0}".format(sys.exc_info()[1])
                )
        else:
            try:
                os.remove(file)  # if database existed but files were deleted
            except PermissionError:
                warn(
                    "Database index could not be updated: {0}".format(sys.exc_info()[1])
                )
            except FileNotFoundError:
                pass

    def find_duplicates(self, columns=None):
        """Find spectra with same conditions. The first duplicated spectrum
        will be ``'False'``, the following will be ``'True'`` (see
        .duplicated()).

        Parameters
        ----------
        columns: list, or ``None``
            columns to find duplicates on. If ``None``, use all conditions.

        Examples
        --------
        >>> db.find_duplicates(columns={'x_e', 'x_N_II'})

            Out[34]:
            file
            20180710_101.spec    True
            20180710_103.spec    True
            dtype: bool

        You can see more examples on the :ref:`Spectrum Database section <label_spectrum_database>`
        of the website.
        """
        dg = self.see(columns=columns).astype(str).duplicated()
        # need to convert eveything as a str to avoid comparaison problems (Minou)
        if columns is None:
            columns = "all"

        if self.verbose:
            print(
                (
                    "{0} duplicate(s) found".format(dg.sum())
                    + " based on columns: {0}".format(columns)
                )
            )

        onlyDuplicatedFiles = dg[dg == True]
        return onlyDuplicatedFiles

    def add(
        self, spectrum: Spectrum, store_name=None, if_exists_then="increment", **kwargs
    ):
        """Add Spectrum to database, whether it's a
        :py:class:`~radis.spectrum.spectrum.Spectrum` object or a file that
        stores one. Check it's not in database already.

        Parameters
        ----------
        spectrum: :py:class:`~radis.spectrum.spectrum.Spectrum` object, or path to a .spec file (``str``)
            if a :py:class:`~radis.spectrum.spectrum.Spectrum` object:
            stores it in the database (using the :py:meth:`~radis.spectrum.spectrum.Spectrum.store`
            method), then adds the file to the database folder.
            if a path to a file (``str``): first copy the file to the database folder,
            then loads the copied file to the database.

        Other Parameters
        ----------------
        store_name: ``str``, or ``None``
            name of the file where the spectrum will be stored. If ``None``,
            name is generated automatically from the Spectrum conditions (see
            ``add_info=`` and ``if_exists_then=``)
        if_exists_then: ``'increment'``, ``'replace'``, ``'error'``, ``'ignore'``
            what to do if file already exists. If ``'increment'`` an incremental digit
            is added. If ``'replace'`` file is replaced (!). If ``'ignore'`` the
            Spectrum is not added to the database and no file is created.
            If ``'error'`` (or anything else) an error is raised. Default ``'increment'``.

        **kwargs: **dict
            extra parameters used in the case where spectrum is a file and a .spec object
            has to be created (useless if `spectrum` is a file already). kwargs
            are forwarded to Spectrum.store() method. See the :meth:`~radis.spectrum.spectrum.Spectrum.store`
            method for more information.

        *Other :py:meth:`~radis.spectrum.spectrum.Spectrum.store` parameters can be given
        as kwargs arguments:*

        compress: 0, 1, 2
            if ``True`` or 1, save the spectrum in a compressed form

            if 2, removes all quantities that can be regenerated with :py:meth:`~radis.spectrum.spectrum.Spectrum.update`,
            e.g, transmittance if abscoeff and path length are given, radiance if
            emisscoeff and abscoeff are given in non-optically thin case, etc.
            If not given, use the value of ``SpecDatabase.binary``
            The performances are usually better if compress = 2. See https://github.com/radis/radis/issues/84.
        add_info: list
            append these parameters and their values if they are in conditions
            example::

                nameafter = ['Tvib', 'Trot']

        discard: list of str
            parameters to exclude. To save some memory for instance
            Default ``['lines', 'populations']``: retrieved Spectrum will loose the
            :meth:`~radis.spectrum.spectrum.Spectrum.line_survey` and
            :meth:`~radis.spectrum.spectrum.Spectrum.plot_populations` methods
            (but it saves a ton of memory!).

        Examples
        --------
        ::

            from radis.tools import SpecDatabase
            db = SpecDatabase(r"path/to/database")     # create or loads database
            db.add(s, discard=['populations'])

        You can see more examples on the :ref:`Spectrum Database section <label_spectrum_database>`
        of the website.

        .. minigallery:: radis.SpecList.add

        See Also
        --------

        :meth:`~radis.tools.database.SpecList.get`,
        :meth:`~radis.tools.database.SpecList.get_unique`,
        :meth:`~radis.tools.database.SpecList.get_closest`,
        :meth:`~spectrum.spectrum.Spectrum.update`
        """

        # Check inputs
        if "path" in kwargs:
            raise ValueError(
                "path is an invalid Parameter. The database path " + "is used"
            )
        compress = kwargs.pop("compress", self.binary)
        store_path = join(self.path, str(store_name))
        if exists(store_path) and if_exists_then == "ignore":
            if self.verbose:
                print(("File already exists : {0}. Ignoring".format(store_name)))
            return

        # First, store the spectrum in a file
        # ... input is a Spectrum. Store it in database and load it from there
        if isinstance(spectrum, Spectrum):
            # add defaults
            if store_name == None:
                if not "add_info" in kwargs:
                    kwargs["add_info"] = self.add_info
                if not "add_date" in kwargs:
                    kwargs["add_date"] = self.add_date
                file = spectrum.store(
                    self.path,
                    compress=compress,
                    if_exists_then=if_exists_then,
                    **kwargs,
                )
            else:
                file = spectrum.store(
                    store_path,
                    compress=compress,
                    if_exists_then=if_exists_then,
                    **kwargs,
                )

        # ... input is a file name. Copy it in database and load it
        elif isinstance(spectrum, str):
            if not exists(spectrum):
                raise FileNotFoundError("File doesnt exist: {0}".format(spectrum))

            fd, name = split(spectrum)

            # Assert a similar case name is not in database already
            if spectrum in list(self.df["file"]):
                raise ValueError(
                    "File already in database: {0}. Database filenames should be unique".format(
                        spectrum
                    )
                )

            if abspath(fd) != abspath(self.path):
                # Assert file doesnt exist in database already
                if name in os.listdir(self.path):
                    raise ValueError(
                        "File already in database folder: {0}".format(name)
                        + ". Use db.update() if you added it there manually"
                    )
                # Ok. Copy it.
                file = join(self.path, name)
                copy2(spectrum, file)

            else:
                warn(
                    "You are manually adding a file that is in the database folder directly. Consider using db.update()"
                )
                file = spectrum

        else:
            raise ValueError("Unvalid Spectrum type: {0}".format(type(spectrum)))

        # Now register the spectrum in the database :
        spectrum_conditions = self._load_new_file(file, binary=compress)
        self.df = self.df.append(spectrum_conditions, ignore_index=True)

        # Update index .csv
        self.print_index()

        return file

    def interpolate(self, **kwconditions):
        """Interpolate existing spectra from the database to generate a new spectrum with conditions kwargs

        Examples
        --------
        ::

            db.interpolate(Tgas=300, mole_fraction=0.3)

        """

        # Not interpolated conditions :
        not_interpolated_conditions = [
            c for c in self.conditions() if c not in kwconditions
        ]
        not_interpolated_conditions = [
            c
            for c in not_interpolated_conditions
            if c not in ["file", "last_modified", "name"]
        ]
        df = self.see(columns=not_interpolated_conditions)

        # assert df values are unique :
        for c in df.columns:
            if df[c].nunique() > 1:
                raise ValueError(
                    f"Spectra in database {self.name} have different values for `{c}`. All conditions except from the one we are interpolating on should be the same"
                )
                # TODO : implement a way to get a subset of a SpecDatabase, so we can do something like db.take(molecule='C2').interpolate(... )

        if len(kwconditions) > 1:
            raise NotImplementedError(
                "Interpolation of multiple conditions is not implemented yet"
            )

        cond = list(kwconditions.keys())[0]
        val = kwconditions[cond]

        cond_list = self.see(columns=[cond]).values[:, 0]
        b = np.argsort(cond_list)

        cond_list = cond_list[b]
        spectra = (self.see().index.values)[b]

        # Interpolate val over cond_list
        pos = np.interp(val, cond_list, np.arange(cond_list.size))
        index = pos.astype(int)
        weight = pos - index

        s_left = self.get_unique(file=spectra[index])
        s_right = self.get_unique(file=spectra[index + 1])

        # linear algebra on spectra directly :
        s_interp = (1 - weight) * s_left + weight * s_right
        # TODO : Will raise an error if multiple values in database. In this
        # case, use s_left.take(var)   with 'var' given in parameters of interpolate()
        # Or tell user to generate a subdatabase with only one spectral array

        s_interp.conditions[
            "interpolated_from"
        ] = f"{spectra[index]}, {spectra[index+1]}"

        return s_interp

    def fit_spectrum(
        self,
        s_exp,
        residual=None,
        normalize=False,
        normalize_how="max",
        conditions="",
        **kwconditions,
    ):
        """Returns the Spectrum in the database that has the lowest residual
        with ``s_exp``.

        Parameters
        ----------
        s_exp: Spectrum
            :class:`~radis.spectrum.spectrum.Spectrum` to fit (typically:
            experimental spectrum)

        Other Parameters
        ----------------
        residual: func, or ``None``
            which residual function to use. If ``None``, use
            :func:`~radis.spectrum.compare.get_residual` with option
            ``ignore_nan=True`` and options ``normalize`` and ``normalize_how``
            as defined by the user.

            ``get_residual`` should have the form::

                lambda s_exp, s, normalize: func(s_exp, s, normalize=normalize)

            where the output is a float.
            Default ``None``
        conditions, **kwconditions: str, **dict
            restrain fitting to only Spectrum that match the given conditions
            in the database. See :meth:`~radis.tools.database.SpecList.get`
            for more information.
        normalize: bool, or Tuple
            see :func:`~radis.spectrum.compare.get_residual`
        normalize_how: 'max', 'area'
            see :func:`~radis.spectrum.compare.get_residual`


        Returns
        -------
        s_best: Spectrum
            closest Spectrum to ``s_exp``

        Examples
        --------
        Using a customized residual function (below: to get the transmittance)::

            from radis import get_residual
            db = SpecDatabase('...')
            db.fit_spectrum(s_exp, get_residual=lambda s_exp, s: get_residual(s_exp, s, var='transmittance'))

        You can see more examples on the :ref:`Spectrum Database section <label_spectrum_database>`
        of the website. More advanced tools for interactive fitting of multi-dimensional, multi-slabs
        spectra can be found in :py:mod:`fitroom`.

        See Also
        --------

        :py:mod:`fitroom`
        """

        if residual is None:
            from radis.spectrum.compare import get_residual

            if len(s_exp.get_vars()) != 1:
                raise ValueError(
                    "Multiple variables in fitted Spectrum. Please "
                    + "define which residual to use with, for instance: "
                    + "`get_residual=lambda s_exp, s: get_residual(s_exp, s, var=SOMETHING)`)"
                )

            def residual(s_exp, s, normalize):
                return get_residual(
                    s_exp,
                    s,
                    var=s_exp.get_vars()[0],
                    ignore_nan=True,
                    normalize=normalize,
                    normalize_how=normalize_how,
                )

        kwconditions.update({"inplace": True})  # dont copy Spectrum to be faster
        spectra = self.get(conditions=conditions, **kwconditions)
        res = np.array([residual(s_exp, s, normalize=normalize) for s in spectra])
        assert not np.isnan(res).any()

        i = np.argmin(np.array(res))

        return spectra[i].copy()  # dont forget to copy the Spectrum we return

    def _load_new_files(self, files, update_register_only=False):
        """Parse files and generate a database.

        Other Parameters
        ----------------
        update_register_only: bool
            if ``True``, load files and update csv but do not keep the Spectrum in memory.
            Default ``False``

        Returns
        -------
        db: pandas.DataFrame
            dataFrame of loaded spectra.

        Notes
        -----
        Can be loaded in parallel using joblib by setting the `nJobs` and `batch_size`
        attributes of :class:`~radis.tools.database.SpecDatabase`.
        See :class:`joblib.parallel.Parallel` for information on the arguments
        """
        db = []

        # get joblib parallel parameters
        nJobs = self.nJobs
        batch_size = self.batch_size
        minimum_nfiles = self.minimum_nfiles

        def funLoad(f):
            return self._load_new_file(
                f, binary=self.binary, update_register_only=update_register_only
            )

        # Sequential loading
        if nJobs == 1 or len(files) < minimum_nfiles:
            if self.verbose:
                print(
                    "*** Loading the database with 1 processor "
                    + "({0} files)***".format(len(files))
                )
            for f in files:
                db.append(funLoad(f))
        # Parallel loading
        else:
            if self.verbose:
                print(
                    "*** Loading the database with {0} ".format(nJobs)
                    + "processor(s) ({0} files)***".format(len(files))
                )

            db = Parallel(
                n_jobs=nJobs, batch_size=batch_size, verbose=5 * self.verbose
            )(delayed(funLoad)(f) for f in files)
        return pd.DataFrame(db)

    def _load_new_file(self, file, binary=False, update_register_only=False):
        """Load spectrum and return Spectrum attributes for insertion in database.

        The Spectrum itself is stored under the "Spectrum" key, and the filename
        under "file".

        Other Parameters
        ----------------
        update_register_only: bool
            if ``True``, load files and update csv but do not keep the Spectrum in memory.
            Default ``False``

        Returns
        -------
        dict: dictionary of conditions of loaded spectrum"""

        s = load_spec(file, binary=binary)
        if self.verbose:
            print(("loaded {0}".format(basename(file))))

        out = s.get_conditions().copy()

        # Add filename, name and a link to the Spectrum object itself
        if s.name == None:
            out["name"] = s.get_name()
        else:
            out["name"] = s.name
        out["file"] = basename(file)
        if update_register_only:
            out["Spectrum"] = None
        else:
            out["Spectrum"] = s
        out["last_modified"] = os.path.getmtime(file)

        return out

    def update_conditions(self):
        """Reloads conditions of all Spectrum in database."""

        # Fetch new conditions, including file and Spectrum object itself
        new_conditions_list = []
        for _, r in self.df.iterrows():
            s = r.Spectrum
            new_conditions = s.get_conditions().copy()
            new_conditions["file"] = r.file
            new_conditions["Spectrum"] = r.Spectrum
            new_conditions_list.append(new_conditions)

        # update DataFrame
        self.df = pd.DataFrame(new_conditions_list)

    def to_dict(self):
        """Returns all Spectra in database under a dictionary, indexed by file.

        Returns
        -------

        out: dict
            {path : Spectrum object} dictionary

        Note
        ----

        ``SpecList.items().values()`` is equivalent to ``SpecList.get()``


        See Also
        --------

        :meth:`~radis.tools.database.SpecList.get`,
        :meth:`~radis.tools.database.SpecList.keys`,
        :meth:`~radis.tools.database.SpecList.values`,
        :meth:`~radis.tools.database.SpecList.items`
        """

        return dict(list(zip(self.df.file, self.df.Spectrum)))


def read_conditions_file(path, verbose=True):
    """Read .csv file with calculation/measurement conditions of all spectra.

    File must have at least the column "file"

    Parameters
    ----------
    path : csv file
        summary of all spectra conditions.

    Returns
    -------
    None.

    """
    # import csv
    # with open(path, 'r') as csvfile:  # Check if comma or semicolon # Note: dialect may be given directly to read_csv ?
    #     dialect=csv.Sniffer().sniff(csvfile.read(), delimiters=',;')

    df = pd.read_csv(path, float_precision="round_trip")
    # thank you https://stackoverflow.com/questions/36909368/precision-lost-while-using-read-csv-in-pandas

    # Force some types
    for boolvar in [
        "self_absorption",
        "db_use_cached",
        "lvl_use_cached",
        "include_neighbouring_lines",
        "export_lines",
        "export_rovib_fractions",
        "load_energies",
        "thermal_equilibrium",
    ]:
        if boolvar in df:
            if df.dtypes[boolvar] == object:
                from radis.misc.basics import str2bool

                if verbose:
                    print(
                        f"Reading {path} : casting column {boolvar} from string type to boolean"
                    )
                df[boolvar] = df[boolvar].map(str2bool)

    # if only "1" in isotopes they are read as numbers and later get() function fails.
    if "isotope" in df:
        df["isotope"] = df["isotope"].astype(str).str.replace(".0", "", regex=True)

    df["Spectrum"] = [None] * len(df["file"])

    return df


def in_database(smatch, db=".", filt=".spec"):
    """Old function."""
    match_cond = smatch.get_conditions()
    for f in [f for f in os.listdir(db) if f.endswith(filt)]:
        fname = join(db, f)
        s = load_spec(fname)
        if s.get_conditions() == match_cond:
            return True
    return False


# %% Test
if __name__ == "__main__":

    import pytest

    pytest.main(["../test/tools/test_database.py", "-s"])
