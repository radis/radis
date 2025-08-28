"""Molecular database (MDB) class

   * MdbExomol is the MDB for ExoMol

Initial code borrowed from the `Exojax <https://github.com/HajimeKawahara/exojax>`__
code (which you should also have a look at!), by @HajimeKawahara, under MIT License.

"""

import os
import pathlib
import urllib.request
import warnings

import numpy as np

from radis.api.dbmanager import DatabaseManager, get_auto_MEMORY_MAPPING_ENGINE
from radis.db.classes import EXOMOL_MOLECULES, EXOMOL_ONLY_ISOTOPES_NAMES

EXOMOL_URL = "https://www.exomol.com/db/"

import bz2
import re
from urllib.request import HTTPError, urlopen

import pandas as pd
from bs4 import BeautifulSoup

from radis.api.hdf5 import vaexsafe_colname


def exact_molname_exomol_to_simple_molname(exact_exomol_molecule_name):
    """convert the exact molname (used in ExoMol) to the simple molname.

    Args:
        exact_exomol_molecule_name: the exact exomol molecule name

    Returns:
        simple molname

    Examples:
        >>> print(exact_molname_exomol_to_simple_molname("12C-1H4"))
        >>> CH4
        >>> print(exact_molname_exomol_to_simple_molname("23Na-16O-1H"))
        >>> NaOH
        >>> print(exact_molname_exomol_to_simple_molname("HeH_p"))
        >>> HeH_p
        >>> print(exact_molname_exomol_to_simple_molname("trans-31P2-1H-2H")) #not working
        >>> Warning: Exact molname  trans-31P2-1H-2H cannot be converted to simple molname
        >>> trans-31P2-1H-2H
        >>> print(exact_molname_exomol_to_simple_molname("1H-2H-16O"))
        >>> H2O
    """

    try:
        molname_simple = _molname_simple_no_exception(exact_exomol_molecule_name)
        return molname_simple
    except:
        print(
            "Warning: Exact molname ",
            exact_exomol_molecule_name,
            "cannot be converted to simple molname",
        )
        return exact_exomol_molecule_name


def _molname_simple_no_exception(exact_exomol_molecule_name):
    t = exact_exomol_molecule_name.split("-")
    molname_simple = ""
    for ele in t:
        alp = "".join(re.findall(r"\D", ele))
        num0 = re.split("[A-Z]", ele)[1]
        if num0.isdigit():
            num = num0
        else:
            num = ""
        molname_simple = molname_simple + alp + num

    # ExoJAX Issue #528
    if molname_simple == "HHO":
        molname_simple = "H2O"
    return molname_simple


def read_def(deff):
    """Exomol IO for a definition file

    Parameters
    ----------
    deff: definition file

    Returns
    -------
    n_Texp : float
        temperature exponent
    alpha_ref : float
        broadening parameter
    molmass : float
        molecular mass
    numinf : List[float]
        limit points (``[w(0), w(1), ..., w(n)]``, n+1 elements) defining
        the spectral ranges appearing in the name of *.trans.bz2 files
        (``["w(0)-w(1)", "w(1)-w(2)", ..., w(n-1)-w(n)]``, n elements)
    numtag : List[str]
        tag for wavelength ranges.

    Note:
        For some molecules, ExoMol provides multiple trans files. numinf and numtag are the ranges and identifiers for the multiple trans files.

    """

    dat = pd.read_csv(deff, sep="#", names=("VAL", "COMMENT"))
    alpha_ref = None
    # texp = None
    n_Texp = None
    ntransf = 1
    maxnu = 0.0
    quantum_labels = []
    unc = False
    for i, com in enumerate(dat["COMMENT"]):
        if "Default value of Lorentzian half-width" in com:
            alpha_ref = float(dat["VAL"][i])
        elif "Default value of temperature exponent" in com:
            n_Texp = float(dat["VAL"][i])
        elif "No. of transition files" in com:
            ntransf = int(dat["VAL"][i])
        elif "Maximum wavenumber (in cm-1)" in com:
            maxnu = float(dat["VAL"][i])
            # maxnu=20000.0
        elif "Isotopologue mass (Da) and (kg)" in com:
            molmass = float(dat["VAL"][i].split()[0])  # in Da (atomic unit)

        elif "Lifetime availability" in com:
            lifetime = int(dat["VAL"][i]) == 1
        elif "Lande g-factor availability" in com:
            lande = int(dat["VAL"][i]) == 1
        elif "Quantum label" in com:
            quantum_labels.append(dat["VAL"][i].strip(" "))
        elif "Uncertainty availability" in com:
            unc = int(dat["VAL"][i]) == 1

    # SOME DEF FILES CONTAINS ERRORS. THESE ARE THE EXCEPTIONS
    if deff.stem == "1H-35Cl__HITRAN-HCl":
        quantum_labels = ["v"]
        # See https://github.com/HajimeKawahara/exojax/issues/330
    if deff.stem == "16O-1H__MoLLIST":
        quantum_labels = ["e/f", "v", "F1/F2", "Es"]

    if deff.stem == "1H-2H-16O__VTT":
        numinf = wavenumber_range_HDO_VTT()
    else:
        numinf = compute_wavenumber_ranges(ntransf, maxnu)

    numtag = wavenumber_tag(numinf)

    output = {
        "n_Texp": n_Texp,
        "alpha_ref": alpha_ref,
        "molmass": molmass,
        "numinf": numinf,
        "numtag": numtag,
        "quantum_labels": quantum_labels,  # array
        "Landé": lande,  # bool
        "lifetime": lifetime,  # bool
        "unc": unc,  # bool uncertainty of line center availability
    }
    return output


def wavenumber_range_HDO_VTT():
    """wave number range for HDO VTT as an exception

    Returns:
        float: numinf
    """
    numinf = np.array(
        [
            0.0,
            250.0,
            500.0,
            750.0,
            1000.0,
            1500.0,
            2000.0,
            2250.0,
            2750.0,
            3500.0,
            4500.0,
            5500.0,
            7000.0,
            9000.0,
            14000.0,
            20000.0,
            26000.0,
        ]
    )
    return numinf


def compute_wavenumber_ranges(ntransf, maxnu):
    """wavenumber range for general case

    Args:
        ntransf: number of the trans files
        maxnu: maximum wavenumber

    Returns:
        float: numinf
    """
    if ntransf > 1:
        dnufile = maxnu / ntransf
        numinf = dnufile * np.array(range(ntransf + 1))
        wavenumber_tag(numinf)
    else:
        numinf = None
    return numinf


def wavenumber_tag(numinf):
    """convert numinf to numtag (used in the name of trans file)


    Args:
        numinf (_type_): wavenumber values for the range

    Returns:
        float: numtag wavenumber tag (such as 00000-00500)
    """
    if numinf is None:
        return ""
    numtag = []
    for i in range(len(numinf) - 1):
        imin = f"{int(numinf[i]):05}"
        imax = f"{int(numinf[i + 1]):05}"
        numtag.append(imin + "-" + imax)
    return numtag


def read_pf(pff):
    """Exomol IO for partition file

    Note:
        T=temperature
        QT=partition function

    Args:
        pff: partition file
    Returns:
        partition data in pandas DataFrame

    """
    dat = pd.read_csv(pff, sep=r"\s+", names=("T", "QT"))
    return dat


def read_trans(transf, engine="vaex"):
    """Exomol IO for a transition file

    Notes
    -----

    Transf format ::

        i_upper=Upper state counting number
        i_lower=Lower state counting number
        A=Einstein coefficient in s-1
        nu_lines=transition wavenumber in cm-1

    See Table 12 in https://arxiv.org/pdf/1603.05890.pdf [Exomol-2016]_

    Parameters
    ----------
    transf: transition file
    engine: parsing engine to use ('vaex', 'csv')

    Returns
    -------
    transition data in vaex/pandas DataFrame

    """
    if engine == "vaex":
        import vaex

        try:  # bz2 compression
            dat = vaex.from_csv(
                transf,
                compression="bz2",
                sep=r"\s+",
                names=("i_upper", "i_lower", "A", "nu_lines"),
                convert=False,  #  file is created by MdbMol
            )
        except:
            try:
                dat = vaex.read_csv(
                    transf,
                    sep=r"\s+",
                    names=("i_upper", "i_lower", "A", "nu_lines"),
                    convert=False,  #  file is created by MdbMol
                )
            except Exception as err:
                raise Exception(
                    f"Error reading {transf}. Either the file does not exist on exomol.com, or was corrupted? You can try to delete it."
                ) from err
    elif engine == "csv":
        try:  # bz2 compression
            dat = pd.read_csv(
                transf,
                compression="bz2",
                sep=r"\s+",
                names=("i_upper", "i_lower", "A", "nu_lines"),
            )
        except:
            try:
                dat = pd.read_csv(
                    transf, sep=r"\s+", names=("i_upper", "i_lower", "A", "nu_lines")
                )
            except Exception as err:
                raise Exception(
                    f"Error reading {transf}. Either the file does not exist on exomol.com, or was corrupted? You can try to delete it."
                ) from err
    else:
        raise NotImplementedError(engine)

    return dat


def read_states(
    statesf, dic_def, engine="vaex", skip_optional_data=True, print_states=False
):
    """Exomol IO for a state file

    Notes
    -----

    States f format ::

        i=state counting number
        E=state energy
        g=state degeneracy
        J=total angular momentum

    See Table 11 in https://arxiv.org/pdf/1603.05890.pdf

    Parameters
    ----------
    statesf: str
        state file
    dic_def: dict
        Info from def file to read extra quantum numbers
    engine: str
        parsing engine to use ('vaex', 'csv')
    skip_optional_data: bool
        If False, fetch all fields which are marked as available in the ExoMol definition
        file. If True, load only the first 4 columns of the states file
        ("i", "E", "g", "J"). The structure of the columns above 5 depend on the
        the definitions file (*.def) and the Exomol version.
        If ``skip_optional_data=False``, two errors may occur:

            - a field is marked as present/absent in the *.def field but is
              absent/present in the *.states file (ie both files are inconsistent).
            - in the updated version of Exomol, new fields have been added in the
              states file of some species. But it has not been done for all species,
              so both structures exist. For instance, the states file of
              https://exomol.com/data/molecules/HCl/1H-35Cl/HITRAN-HCl/ follows the
              structure described in [1]_, unlike the states file of
              https://exomol.com/data/molecules/NO/14N-16O/XABC/ which follows the
              structure described in [2]_.
    print_states: bool
        print some info. When .def file looks inconsistent with .states file, turn ON and report them in Issue.


    Returns
    -------
    states data in pandas DataFrame

    If ``'vaex'``, also writes a local hdf5 file `statesf.with_suffix('.hdf5')`


    References
    ----------

    .. [1] Tennyson, J., Yurchenko, S. N., Al-Refaie, A. F., Barton, E. J., Chubb, K. L., Coles, P. A., … Zak, E. (2016). The ExoMol database: molecular line lists for exoplanet and other hot atmospheres. https://doi.org/10.1016/j.jms.2016.05.002
    .. [2] Tennyson, J., Yurchenko, S. N., Al-Refaie, A. F., Clark, V. H. J., Chubb, K. L., Conway, E. K., … Yurchenko, O. P. (2020). The 2020 release of the ExoMol database: Molecular line lists for exoplanet and other hot atmospheres. Journal of Quantitative Spectroscopy and Radiative Transfer, 255, 107228. https://doi.org/10.1016/j.jqsrt.2020.107228

    """
    # we read first 4 columns for ("i", "E", "g", "J"),
    # skip lifetime, skip Landé g-factor,
    # read quantum numbers
    N_otherfields = np.size(dic_def["quantum_labels"])
    quantum_labels = dic_def["quantum_labels"]

    mandatory_fields = ("i", "E", "g", "J")
    N_mandatory_fields = len(mandatory_fields)
    mandatory_usecol = np.arange(N_mandatory_fields)
    if skip_optional_data:
        usecol = mandatory_usecol
        names = mandatory_fields
    else:
        label = np.array(["unc", "lifetime", "Lande"])
        mask = np.array([dic_def["unc"], dic_def["Landé"], dic_def["lifetime"]])
        N_except = np.sum(mask)
        usecol = np.arange(0, N_mandatory_fields + N_except + N_otherfields)
        names = mandatory_fields + tuple(label[mask]) + tuple(quantum_labels)
        # The definitions file (*.def) specifies which fields are available
        # in the states fiel (*.states). Check the number of available fields
        # (see *.def) match the numbers of columns in the states file (*.states).
        # Otherwise, both files are inconsistent.
        with bz2.open(statesf, "rb") as f:
            firstline = f.readline().decode("utf-8")  # read 1dt line
            splitline = [x for x in re.split(r"\s+", firstline) if x != ""]
            # splitline = re.split(r"\s+", firstline)
            f.close()
        if len(splitline) != len(names):
            warnings.warn(
                "There appear to be further additional columns in .state file,"
                + "which are not defined in .def file. We skip them."
            )

    if engine == "vaex":
        import vaex

        # TODO Refactor: move in DataFileManager
        try:
            dat = vaex.from_csv(
                statesf,
                compression="bz2",
                sep=r"\s+",
                usecols=usecol,
                names=names,
                convert=False,  # written in MolDB
            )
        except:
            try:
                dat = vaex.read_csv(
                    statesf,
                    sep=r"\s+",
                    usecols=usecol,
                    names=names,
                    convert=False,  # written in MolDB
                )

            except Exception as err:
                raise Exception(
                    f"Error reading {statesf}. Either the file does not exist on exomol.com, or was corrupted? You can try to delete it."
                ) from err
    elif engine == "csv":
        try:
            dat = pd.read_csv(
                statesf, compression="bz2", sep=r"\s+", usecols=usecol, names=names
            )
        except:  #!!!TODO What was the expected error?
            dat = pd.read_csv(statesf, sep=r"\s+", usecols=usecol, names=names)
    else:
        raise NotImplementedError(engine)
    if print_states:
        print(dat)
    return dat


def pickup_gE(states, trans, dic_def, skip_optional_data=True, engine="vaex"):
    """extract g_upper (gup), E_lower (elower), and J_lower and J_upper from states
    DataFrame and insert them into the transition DataFrame.

    Parameters
    ----------
    states: states DataFrame  - the i, E, g, J are in the 4 first columns
    trans: transition numpy array
    dic_def: Informations about additional quantum labels
    skip_optional_data: bool . If True fetch all quantum labels in dic_def['quantum_labels'] from states into transitions (_l for lower, _u for upper states)

    Returns
    -------
    A, nu_lines, elower, gup, jlower, jupper, mask, **quantum_labels

    Notes
    -----
    We first convert pandas DataFrame to ndarray. The state counting numbers in states DataFrame is used as indices of the new array for states (newstates). We remove the state count numbers as the column of newstate, i.e. newstates[:,k] k=0: E, 1: g, 2: J. Then, we can directly use the state counting numbers as mask.

    States by default has columns ::

        #       i      E             g    J    v
        0       1      0.0           1    0    0
        1       2      1.448467      3    1    0
        2       3      4.345384      5    2    0
        3       4      8.690712      7    3    0


    """

    ### Step 1. Essential quantum number for spectra
    # ----------------------------------------------

    def map_add(trans, col, new_col, trans_key, states_key="i"):
        """Lookup `key` in states and add it in trans, using the level ``i``
        for upper and lower state

        Examples
        --------
        ::

            trans = map_add(trans, "E", "E_lower", "i_lower")
        """
        if engine == "pytables":
            # Rename the columns in the states DataFrame
            states_map = states.copy()
            states_map.rename(
                columns={col: new_col, states_key: trans_key}, inplace=True
            )
            states_map = states_map[[new_col, trans_key]]  # drop useless columns
            trans = trans.join(states_map.set_index(trans_key), on=trans_key)
        elif engine == "vaex":
            col = vaexsafe_colname(col)
            new_col = vaexsafe_colname(new_col)

            # WIP. First implementation with join(). Use Map() will probably be faster.
            trans.join(
                states[states_key, col],
                left_on=trans_key,
                right_on=states_key,
                inplace=True,
            )
            trans.drop(states_key, inplace=True)
            trans.rename(col, new_col)
        return trans

    trans = map_add(trans, "g", "gup", "i_upper")
    trans = map_add(trans, "J", "jlower", "i_lower")
    trans = map_add(trans, "J", "jupper", "i_upper")
    trans = map_add(trans, "E", "elower", "i_lower")

    def has_nan(column):
        try:  # Vaex
            return column.countnan() > 0
        except AttributeError:  # Pandas
            return column.hasnans

    if not "nu_lines" in trans or has_nan(trans["nu_lines"]):
        trans = map_add(trans, "E", "eupper", "i_upper")
        trans["nu_lines"] = trans["eupper"] - trans["elower"]

    ### Step 2. Extra quantum numbers (e/f parity, vib and rot numbers)
    # -----------------------------------------------------------------
    if not skip_optional_data:
        for q in dic_def["quantum_labels"]:
            trans = map_add(trans, q, f"{q}_l", "i_lower")
            trans = map_add(trans, q, f"{q}_u", "i_upper")

    return trans


# def pickup_gEslow(states, trans):
#     """Slow version to extract g_upper (gup) and E_lower (elower) from states DataFrame and insert them to transition DataFrame.

#     Note:
#        This code is the (thousands times) slower version of pickup_gE. However, we did not remove this one just because this version is much easier to understand.

#     Args:
#        states: states pandas DataFrame
#        trans: transition pandas DataFrame

#     Returns:
#        A, nu_lines, elower, gup

#     """
#     import tqdm

#     E = states["E"].values
#     g = states["g"].values

#     # insert new columns in transition array
#     trans["gup"] = 0
#     trans["elower"] = 0.0

#     for k, i in tqdm.tqdm(enumerate(states["i"])):
#         # transition upper state
#         mask_upper = trans["i_upper"] == i
#         trans["gup"][mask_upper] = g[k]
#         # transition lower state
#         mask_lower = trans["i_lower"] == i
#         trans["elower"][mask_lower] = E[k]

#     A = trans["A"].to_numpy()
#     nu_lines = trans["nu_lines"].to_numpy()
#     elower = trans["elower"].to_numpy()
#     gup = trans["gup"].to_numpy()
#     return A, nu_lines, elower, gup


def read_broad(broadf, output="pytables"):
    """Reading broadening file (.broad)
    Parameters
    ----------
    broadf: .broad file
    engine: "pytables" (default) for NumPy/Pandas processing, "vaex" for Vaex processing.

    Returns
    -------
    broadening info in bdat form (pandas or vaex dataframe), defined by this instance.

    Notes
    -----
    See Table 16 in https://arxiv.org/pdf/1603.05890.pdf
    """
    column_names = [
        "code",
        "alpha_ref",
        "n_Texp",
        "jlower",
        "jupper",
        "kalower",
        "kclower",
        "kaupper",
        "kcupper",
        "v1lower",
        "v2lower",
        "v3lower",
        "v1upper",
        "v2upper",
        "v3upper",
    ]

    if output == "vaex":
        import vaex

        bdat = vaex.from_csv(broadf, sep=r"\s+", names=column_names)
    else:
        bdat = pd.read_csv(broadf, sep=r"\s+", names=column_names)

    return bdat


def check_code_level(bdat, output="pytables"):
    """checking code level in .broad
    Args:
        bdat: exomol .broad data given by exomolapi.read_broad
        output: "pytables" (default) for NumPy/Pandas processing, "vaex" for Vaex processing.
    Returns:
        code level: None, a0, a1, other codes unavailable currently,
        if a0 and a1 are available, a1 is returned.
        the other cases returns None
    """
    if output == "vaex":
        input_array = bdat.unique(bdat.code)
    else:
        if type(bdat) == dict:
            input_array = np.unique(bdat["code"])
        else:
            input_array = pd.unique(bdat.code)
    if np.array_equal(input_array, np.array(["a0"])):
        return "a0"
    elif np.array_equal(input_array, np.array(["a1"])):
        return "a1"
    elif np.array_equal(np.sort(input_array), np.array(["a0", "a1"])):
        return "a1"
    elif np.array_equal(input_array, np.array(["m0"])):
        return "m0"
    return None


def make_j2b(
    bdat, alpha_ref_default=0.07, n_Texp_default=0.5, jlower_max=None, output="pytables"
):
    """compute j2b (code a0, map from jlower to alpha_ref)

    Args:
        bdat: exomol .broad data given by exomolapi.read_broad
        alpha_ref_default: default value
        n_Texp_default: default value
        jlower_max: maximum number of jlower
        engine: "pytables" (default) for NumPy/Pandas processing, "vaex" for Vaex processing.
    Returns:
        j2alpha_ref[jlower] provides alpha_ref for jlower
        j2n_Texp[jlower]  provides nT_exp for jlower
    """
    # a0
    if output == "vaex":
        import vaex

        bdat = bdat[bdat["code"] == "a0"]
        jlower_arr = bdat["jlower"].values
        alpha_ref_arr = bdat["alpha_ref"].values
        n_Texp_arr = bdat["n_Texp"].values
    else:
        cmask = bdat["code"] == "a0"
        jlower_arr = np.array(bdat["jlower"][cmask], dtype=int)
        alpha_ref_arr = np.array(bdat["alpha_ref"][cmask])
        n_Texp_arr = np.array(bdat["n_Texp"][cmask])

    # Determine the array size based on jlower_max
    if jlower_max is None:
        Nblower = np.nanmax(jlower_arr) + 1
    else:
        Nblower = np.nanmax([jlower_max, np.max(jlower_arr)]) + 1

    if output == "vaex":
        # Initialize arrays with default alpha_ref and n_Texp
        df_defaults = vaex.from_arrays(
            jlower=vaex.vrange(0, Nblower, dtype="int64"),
            alpha_ref=np.full(Nblower, alpha_ref_default),
            n_Texp=np.full(Nblower, n_Texp_default),
        )
        df_defaults = df_defaults.join(bdat, on="jlower", how="left", rsuffix="_new")

        df_defaults["alpha_ref"] = df_defaults["alpha_ref_new"].fillna(
            alpha_ref_default
        )
        df_defaults["n_Texp"] = df_defaults["n_Texp_new"].fillna(n_Texp_default)

        # Populate the mapping arrays using known broadening coefficients
        j2alpha_ref = df_defaults["alpha_ref"].values
        j2n_Texp = df_defaults["n_Texp"].values
    else:
        # Initialize arrays with default alpha_ref and n_Texp
        j2alpha_ref = np.full(int(Nblower), alpha_ref_default)
        j2n_Texp = np.full(int(Nblower), n_Texp_default)

        # Populate the mapping arrays using known broadening coefficients
        j2alpha_ref[jlower_arr] = alpha_ref_arr
        j2n_Texp[jlower_arr] = n_Texp_arr

    # Raise a minor warning if default values are used for high J values
    if Nblower > (np.max(jlower_arr) + 1):
        import warnings

        from radis.misc.warning import AccuracyWarning

        warnings.warn(
            AccuracyWarning(
                f"The default broadening parameter (alpha = {alpha_ref_default} cm^-1 and n = {n_Texp_default}) are used for J'' > {np.max(jlower_arr)} up to J'' = {Nblower}"
            )
        )

    return j2alpha_ref, j2n_Texp


def make_jj2b(bdat, j2alpha_ref_def, j2n_Texp_def, jupper_max=None, output="pytables"):
    """compute jj2b (code a1, map from (jlower, jupper) to alpha_ref and n_Texp)

    Args:
        bdat: exomol .broad data given by exomolapi.read_broad
        j2alpha_ref_def: default value from a0
        j2n_Texp_def: default value from a0
        jupper_max: maximum number of jupper
        engine: "pytables" (default) for NumPy/Pandas processing, "vaex" for Vaex processing.
    Returns:
        jj2alpha_ref[jlower,jupper] provides alpha_ref for (jlower, jupper)
        jj2n_Texp[jlower,jupper]  provides nT_exp for (jlower, jupper)
    Note:
        The pair of (jlower, jupper) for which broadening parameters are not given, jj2XXX contains None.
    """
    # a1
    if output == "vaex":
        import vaex

        bdat = bdat[bdat["code"] == "a1"]
        jlower_arr = bdat["jlower"].values
        jupper_arr = bdat["jupper"].values
        alpha_ref_arr = bdat["alpha_ref"].values
        n_Texp_arr = bdat["n_Texp"].values
    else:
        cmask = bdat["code"] == "a1"
        jlower_arr = np.array(bdat["jlower"][cmask], dtype=int)
        jupper_arr = np.array(bdat["jupper"][cmask], dtype=int)
        alpha_ref_arr = np.array(bdat["alpha_ref"][cmask])
        n_Texp_arr = np.array(bdat["n_Texp"][cmask])

    # Determine the array size based on jupper_max
    if jupper_max is None:
        Nbupper = np.max(jupper_arr) + 1
    else:
        Nbupper = np.max([jupper_max, np.max(jupper_arr)]) + 1

    if output == "vaex":
        df = vaex.from_arrays(
            j2alpha_ref_def=[j2alpha_ref_def], j2n_Texp_def=[j2n_Texp_def]
        )
        df["j2alpha_ref_column"] = df["j2alpha_ref_def"]
        df["j2n_Texp_column"] = df["j2n_Texp_def"]

        df["jj2alpha_ref"] = df["j2alpha_ref_column"] * 1
        df["jj2n_Texp"] = df["j2n_Texp_column"] * 1
        ones_array = np.ones(Nbupper)

        jj2alpha_ref = np.outer(df["j2alpha_ref_def"].to_numpy(), ones_array)
        jj2n_Texp = np.outer(df["j2n_Texp_def"].to_numpy(), ones_array)
    else:
        jj2alpha_ref = j2alpha_ref_def[:, np.newaxis] * np.ones(Nbupper)
        jj2n_Texp = j2n_Texp_def[:, np.newaxis] * np.ones(Nbupper)

    jj2alpha_ref[jlower_arr, jupper_arr] = alpha_ref_arr
    jj2n_Texp[jlower_arr, jupper_arr] = n_Texp_arr

    return jj2alpha_ref, jj2n_Texp


def make_j2b_m0(bdat, alpha_ref_default=0.07, n_Texp_default=0.5, jlower_max=None):
    """compute j2b (code m0, map from |m| to alpha_ref)

    Args:
        bdat: exomol .broad data given by exomolapi.read_broad
        alpha_ref_default: default value
        n_Texp_default: default value
        jlower_max: maximum number of jlower
    Returns:
        j2alpha_ref[jlower] provides alpha_ref for jlower
        j2n_Texp[jlower]  provides nT_exp for jlower
    """
    # m0
    cmask = bdat["code"] == "m0"
    jlower_arr = np.array(bdat["jlower"][cmask], dtype=int)
    alpha_ref_arr = np.array(bdat["alpha_ref"][cmask])
    n_Texp_arr = np.array(bdat["n_Texp"][cmask])

    # Determine the array size based on jlower_max
    if jlower_max is None:
        Nblower = np.max(jlower_arr) + 1
    else:
        Nblower = np.max([jlower_max, np.max(jlower_arr)]) + 1

    # Initialize arrays with default alpha_ref and n_Texp
    j2alpha_ref = np.full(Nblower, alpha_ref_default)
    j2n_Texp = np.full(Nblower, n_Texp_default)

    # Populate the mapping arrays using known broadening coefficients
    j2alpha_ref[jlower_arr] = alpha_ref_arr
    j2n_Texp[jlower_arr] = n_Texp_arr

    # Raise a minor warning if default values are used for high J values
    if Nblower > (np.max(jlower_arr) + 1):
        import warnings

        from radis.misc.warning import AccuracyWarning

        warnings.warn(
            AccuracyWarning(
                f"The default broadening parameter (alpha = {alpha_ref_default} cm^-1 and n = {n_Texp_default}) are used for J'' > {np.max(jlower_arr)} up to J'' = {Nblower}"
            )
        )

    return j2alpha_ref, j2n_Texp


def get_exomol_full_isotope_name(molecule, isotope):
    """Get full isotope name for ``molecule`` and ``isotope`` number

    Parameters
    ----------
    molecule: str
    isotope: int
        terrestrial abundance

    Examples
    --------
    ::

        get_exomol_full_isotope_name("CH4", 1)
        >>> '12C-1H4'

    See Also
    --------
    :py:func:`~radis.api.exomolapi.get_exomol_database_list`"""

    if molecule not in EXOMOL_MOLECULES:
        raise ValueError(
            f"Molecule {molecule} not available in ExoMol. Change database, or choose one of ExoMol available molecules: {sorted(EXOMOL_MOLECULES)}"
        )

    if (molecule, isotope) in EXOMOL_ONLY_ISOTOPES_NAMES:
        return EXOMOL_ONLY_ISOTOPES_NAMES[(molecule, isotope)]
    else:
        # Read and convert from HITRAN molecules
        from radis.db.molparam import MolParams

        mp = MolParams()
        return mp.get(molecule, isotope, "isotope_name_exomol")


def get_list_of_known_isotopes(molecule):
    """find all isotopes until error ensues"""

    i = 1
    isotope_list = []
    while True:
        try:
            iso_name = get_exomol_full_isotope_name(molecule, i)
        except:
            break
        else:
            isotope_list.append(iso_name)
        finally:
            i += 1
    return isotope_list


def get_exomol_database_list(molecule, isotope_full_name=None):
    """Parse ExoMol website and return list of available databases, and recommended database

    Parameters
    ----------
    molecule: str
    isotope_full_name: str
        isotope full name (ex. ``12C-1H4`` for CH4,1). Get it from
        :py:func:`radis.api.exomolapi.get_exomol_full_isotope_name`

    Returns
    -------
    list of databases, database recommended by ExoMol

    Examples
    --------
    Get CH4 from ExoMol :
    ::

        databases, recommended = get_exomol_database_list("CH4", "12C-1H4")
        >>> ['xsec-YT10to10', 'YT10to10', 'YT34to10'], 'YT34to10'

    Or combine with :py:func:`~radis.api.exomolapi.get_exomol_full_isotope_name` to
    get the isopologue (sorted by terrestrial abundance) ::

        from radis.api.exomolapi import get_exomol_database_list, get_exomol_full_isotope_name
        databases, recommended = get_exomol_database_list("CH4", get_exomol_full_isotope_name("CH4", 1))
        >>> ['xsec-YT10to10', 'YT10to10', 'YT34to10'], 'YT34to10'


    .. minigallery:: radis.api.exomolapi.get_exomol_database_list


    See Also
    --------
    :py:func:`~radis.api.exomolapi.get_exomol_full_isotope_name`
    """

    if isotope_full_name is None:
        raise ValueError(
            f"Give isotope name. List of known isotopes for {molecule} : {get_list_of_known_isotopes(molecule)}"
        )

    url = f"https://exomol.com/data/molecules/{molecule}/{isotope_full_name}"
    try:
        response = urlopen(url).read()
    except HTTPError as err:
        if isotope_full_name not in get_list_of_known_isotopes(molecule):
            extra = f". Isotope name {isotope_full_name} is not in list of known isotopes : {get_list_of_known_isotopes(molecule)}"
        else:
            extra = ""
        raise ValueError(f"HTTPError opening url={url}" + extra) from err

    soup = BeautifulSoup(
        response, features="lxml"
    )  # make soup that is parse-able by bs

    # Recommended database
    rows = soup.find_all(
        "a", {"class": "list-group-item link-list-group-item recommended"}
    )
    databases_recommended = [r.get_attribute_list("title")[0] for r in rows]
    databases_recommended = list(np.unique(databases_recommended))

    # All others
    rows = soup.find_all("a", {"class": "list-group-item link-list-group-item"})
    databases = [r.get_attribute_list("title")[0] for r in rows]

    if len(databases_recommended) > 1:
        old_recommended = databases_recommended.copy()

        ### Known exception: SiO cross-section
        # this is a cross-section dataset, shouldn't be used. Reverse and use the other one:
        exception_SiO = (
            isotope_full_name == "28Si-16O"
            and databases_recommended[0] == "xsec-SiOUVenIR"
        )
        if exception_SiO:
            databases_recommended = databases_recommended[::-1]

        ### Known exception: DTU cross-section
        exception_DTU = "DTU" in databases_recommended
        if exception_DTU:
            databases_recommended.remove("DTU")

        if exception_SiO or exception_DTU:
            print(
                f"Multiple recommended databases found for {molecule} in ExoMol : {old_recommended}. {isotope_full_name} is an exception, using {databases_recommended}"
            )
        else:
            databases_recommended[::-1]
            print(
                f"Multiple recommended databases found for {molecule} in ExoMol : {old_recommended}. This is unexpected. Trying with the first: {databases_recommended}"
            )

    databases = databases + databases_recommended

    if len(databases_recommended) > 0:
        recommended_database = databases_recommended[0]
    else:
        recommended_database = False

    return list(np.unique(databases)), recommended_database


# def fetch_exomol_molecule_list():
#     """Parse ExoMol website and return list of available databases, and recommended database

#     Parameters
#     ----------
#     molecule: str
#     isotope_full_name: str
#         isotope full name (ex. ``12C-1H4`` for CH4,1). Get it from
#         :py:func:`radis.api.exomolapi.get_exomol_full_isotope_name`

#     Returns
#     -------

#     Examples
#     --------
#     Get CH4 from ExoMol :
#     ::
#         databases, recommended = get_exomol_database_list("CH4", "12C-1H4")
#         >>> ['xsec-YT10to10', 'YT10to10', 'YT34to10'], 'YT34to10'

#     Or combine with :py:func:`~radis.api.exomolapi.get_exomol_full_isotope_name` to
#     get the isopologue (sorted by terrestrial abundance) ::

#         from radis.api.exomolapi import get_exomol_database_list, get_exomol_full_isotope_name
#         databases, recommended = get_exomol_database_list("CH4", get_exomol_full_isotope_name("CH4", 1))
#         >>> ['xsec-YT10to10', 'YT10to10', 'YT34to10'], 'YT34to10'

#     .. minigallery:: radis.api.exomolapi.get_exomol_database_list

#     See Also
#     --------
#     :py:func:`~radis.api.exomolapi.get_exomol_full_isotope_name`
#     """

# url = f"https://exomol.com/data/molecules/"
# try:
#     response = urlopen(url).read()
# except HTTPError as err:
#     raise ValueError(f"HTTPError opening url={url}") from err

# soup = BeautifulSoup(
#     response, features="lxml"
# )  # make soup that is parse-able by bs

# # Recommended database
# rows = soup.find_all(
#     "a", {"class": "list-group-item link-list-group-item recommended"}
# )
# databases_recommended = [r.get_attribute_list("title")[0] for r in rows]

# # All others
# rows = soup.find_all("a", {"class": "list-group-item link-list-group-item"})
# databases = [r.get_attribute_list("title")[0] for r in rows]

# if len(databases_recommended) > 1:
#     # Known exceptions :
#     if (
#         isotope_full_name == "28Si-16O"
#         and databases_recommended[0] == "xsec-SiOUVenIR"
#     ):
#         # this is a cross-section dataset, shouldn't be used. Reverse and use the other one:
#         databases_recommended = databases_recommended[::-1]
#     else:
#         print(
#             f"Multiple recommended databases found for {molecule} in ExoMol : {databases_recommended}. This is unexpected. Using the first"
#         )

# databases = databases + databases_recommended

# return databases, databases_recommended[0]


class MdbExomol(DatabaseManager):
    """molecular database of ExoMol

    MdbExomol is a class for ExoMol.

    Parameters
    ----------
    path: str
        path for Exomol data directory/tag. For instance, "/home/CO/12C-16O/Li2015"
    nurange: array
        wavenumber range list (cm-1) or wavenumber array
    margin: float
        margin for nurange (cm-1)
    crit: float
        line strength lower limit for extraction
    bkgdatm: str
        background atmosphere for broadening. e.g. H2, He,
    broadf: bool
        if False, the default broadening parameters in .def file is used, default is True
    broadf_download: bool
        if False, not try to download potential broadening files, default is True

    Other Parameters
    ----------------
    engine : str
        which memory mapping engine to use : 'vaex', 'pytables' (HDF5), 'feather'
    skip_optional_data : bool
        If False, fetch all fields which are marked as available in the ExoMol definition
        file. If True, load only the first 4 columns of the states file
        ("i", "E", "g", "J"). The structure of the columns above 5 depend on the
        the definitions file (*.def) and the Exomol version.
        If ``skip_optional_data=False``, two errors may occur:

            - a field is marked as present/absent in the *.def field but is
              absent/present in the *.states file (ie both files are inconsistent).
            - in the updated version of Exomol, new fields have been added in the
              states file of some species. But it has not been done for all species,
              so both structures exist. For instance, the states file of
              https://exomol.com/data/molecules/HCl/1H-35Cl/HITRAN-HCl/ follows the
              structure described in [1]_, unlike the states file of
              https://exomol.com/data/molecules/NO/14N-16O/XABC/ which follows the
              structure described in [2]_.

    Notes
    -----

    The trans/states files can be very large. For the first time to read it,
    we convert it to the feather or hdf5-format. After the second-time,
    we use the feather/hdf5 format instead.

    Examples
    --------
    ::

        # Init database, download files if needed.
        mdb = MdbExomol(
            local_path,
            molecule=molecule,
            name=databank_name,
            local_databases=local_databases,
            # nurange=[load_wavenum_min, load_wavenum_max],
            engine="vaex",
        )

        # Get cache files to load :
        mgr = mdb.get_datafile_manager()
        local_files = [mgr.cache_file(f) for f in mdb.trans_file]

        # Load files
        df = mdb.load(
            local_files,
            columns=columns_exomol,
            lower_bound=([('nu_lines', load_wavenum_min)] if load_wavenum_min else []) + ([("Sij0", mdb.crit)] if not np.isneginf(mdb.crit) else []),
            upper_bound=([('nu_lines', load_wavenum_max)] if load_wavenum_max else []),
            output="jax", # or "pytables", "vaex"
        )

    .. minigallery:: radis.fetch_exomol


    DataFrame columns
    -----------------

    nu_lines (nd array): line center (cm-1)
    Sij0 (nd array): line strength at T=Tref (cm)
    dev_nu_lines (np array): line center in device (cm-1)
    logsij0 (np array): log line strength at T=Tref
    A (np array): Einstein A coefficient
    elower (np array): the lower state energy (cm-1)
    gpp (np array): statistical weight
    jlower (np array): J_lower
    jupper (np array): J_upper
    n_Tref (np array): temperature exponent
    alpha_ref (np array): alpha_ref (gamma0)
    n_Tref_def: default temperature exponent in .def file, used for jlower not given in .broad
    alpha_ref_def: default alpha_ref (gamma0) in .def file, used for jlower not given in .broad

    References
    ----------

    .. [1] Tennyson, J., Yurchenko, S. N., Al-Refaie, A. F., Barton, E. J., Chubb, K. L., Coles, P. A., … Zak, E. (2016). The ExoMol database: molecular line lists for exoplanet and other hot atmospheres. https://doi.org/10.1016/j.jms.2016.05.002
    .. [2] Tennyson, J., Yurchenko, S. N., Al-Refaie, A. F., Clark, V. H. J., Chubb, K. L., Conway, E. K., … Yurchenko, O. P. (2020). The 2020 release of the ExoMol database: Molecular line lists for exoplanet and other hot atmospheres. Journal of Quantitative Spectroscopy and Radiative Transfer, 255, 107228. https://doi.org/10.1016/j.jqsrt.2020.107228

    See also
    --------

    MdbExomol is compatible with Exojax :py:class:`exojax.spec.api.MdbExomol`

    """

    # TODO : inherit from DatabaseManager or similar

    # @dev: In exojax this class is defined in exojax/spec/moldb.py
    # see https://github.com/HajimeKawahara/exojax/blob/develop/src/exojax/spec/moldb.py
    # It uses Jax arrays (jnp). Here RADIS uses Numpy arrays.

    def __init__(
        self,
        path,
        molecule,
        database=None,
        local_databases=None,
        name="EXOMOL-{molecule}",
        nurange=[0.0, np.inf],
        margin=0.0,
        crit=-np.inf,
        bkgdatm="Air",  # TODO: use Air whenever possible (consistent with HITRAN/HITEMP). This is not a parameter for the moment.
        broadf=True,
        broadf_download=True,
        engine="vaex",
        verbose=True,
        cache=True,
        skip_optional_data=True,
    ):
        super().__init__(
            name,
            molecule,
            local_databases,
            engine,
            verbose=verbose,
        )
        assert cache  # cache only used for cache='regen' or cache='force' modes, cache=False is not expected

        if engine == "default":
            from radis import config

            engine = config["MEMORY_MAPPING_ENGINE"]
            if engine == "auto":
                engine = get_auto_MEMORY_MAPPING_ENGINE()
        self.engine = engine

        self.path = pathlib.Path(path)
        if local_databases is not None:
            self.path = pathlib.Path(local_databases).expanduser() / self.path

        t0 = self.path.parents[0].stem
        molec = t0 + "__" + str(self.path.stem)

        self.crit = crit
        self.margin = margin
        self.nurange = [np.min(nurange), np.max(nurange)]
        self.wmin, self.wmax = np.min(nurange), np.max(nurange)
        self.broadf = broadf
        self.broadf_download = broadf_download

        # Where exomol files are
        self.states_file = self.path / pathlib.Path(molec + ".states.bz2")
        self.pf_file = self.path / pathlib.Path(molec + ".pf")
        self.def_file = self.path / pathlib.Path(molec + ".def")

        self.broad_partners = [
            "H2",
            "He",
            "air",
            "self",
            "Ar",
            "CH4",
            "CO",
            "CO2",
            "H2",
            "H2O",
            "N2",
            "NH3",
            "NO",
            "O2",
            "NH3",
            "CS",
        ]  # see ExoMol paper 2024
        self.broad_files = {
            partner: self.path / pathlib.Path(t0 + "__" + partner + ".broad")
            for partner in self.broad_partners
        }

        mgr = self.get_datafile_manager()
        if not self.def_file.exists():
            self.download(molec, extension=[".def"])
        if not self.pf_file.exists():
            self.download(molec, extension=[".pf"])
        if (
            not self.states_file.exists()
            and not mgr.cache_file(self.states_file).exists()
        ):
            self.download(molec, extension=[".states.bz2"])
        # will attempt a download as long as "air" is not present in the database
        if (
            (not self.broad_files["air"].exists())
            and self.broadf
            and self.broadf_download
        ):
            self.download(molec, extension=[".broad"])

        # Add molecule name
        tag = molec.split("__")
        self.isotope_fullname = tag[0]
        self.molecule = exact_molname_exomol_to_simple_molname(tag[0])
        # self.isotope = 1  # Placeholder. TODO : implement parsing of other isotopes.

        # load def
        dic_def = read_def(self.def_file)  # approx. 3 ms
        self.n_Texp_def = dic_def["n_Texp"]
        self.alpha_ref_def = dic_def["alpha_ref"]
        self.molmass = dic_def["molmass"]

        #  default n_Texp value if not given
        if self.n_Texp_def is None:
            self.n_Texp_def = 0.5
            warnings.warn(
                Warning(
                    f"""
                    No default broadening exponent in def file. Assigned n = {self.n_Texp_def}
                    """
                )
            )
        #  default alpha_ref value if not given
        if self.alpha_ref_def is None:
            self.alpha_ref_def = 0.07
            warnings.warn(
                Warning(
                    f"""
                    No default broadening in def file. Assigned alpha_ref = {self.alpha_ref_def}
                    """
                )
            )

        # load states
        if cache == "regen" and mgr.cache_file(self.states_file).exists():
            if self.verbose:
                print("Removing existing file ", mgr.cache_file(self.states_file))
            os.remove(mgr.cache_file(self.states_file))
        if mgr.cache_file(self.states_file).exists():
            states = mgr.read(mgr.cache_file(self.states_file))
        else:
            if cache == "force":
                raise ValueError(
                    f"Cache file {str(mgr.cache_file(self.states_file))} does not exist"
                )
            print(
                f"Note: Caching states data to the {engine} format. After the second time, it will become much faster."
            )
            states = read_states(
                self.states_file,
                dic_def,
                engine="vaex" if engine == "vaex" else "csv",
                skip_optional_data=skip_optional_data,
            )
            mgr.write(mgr.cache_file(self.states_file), states)

        # load pf
        pf = read_pf(self.pf_file)
        self.gQT = pf["QT"].to_numpy()  # grid QT
        self.T_gQT = pf["T"].to_numpy()  # T forgrid QT

        # trans file(s)
        # Compute linestrengths or retrieve them from cache
        self.Tref = 296.0
        self.QTref = np.array(self.QT_interp(self.Tref))

        # Download files

        # Generate list of files
        # ---------------------
        # Case1 : transitions are stored in a single file
        if dic_def["numinf"] is None:
            self.trans_file = [self.path / pathlib.Path(molec + ".trans.bz2")]
            self.num_tag = [None]

        # Case2 : Transitions are stored in multiple files:
        else:  # dic_def["numinf"] is not None
            # dic_def["numinf"] contains the limit points ``[w(0), w(1), ..., w(n)]``
            # (n+1 elements) defining the spectral ranges appearing in the list
            # dic_def["numtag"] which looks like ``["w(0)-w(1)", "w(1)-w(2)", ..., w(n-1)-w(n)]``
            # (n elements)
            # imin :
            # index i of the "w(i)-w(i+1)" element in dic_def["numtag"] such
            # that nurange[0]<=w(i)
            imin = (
                np.searchsorted(dic_def["numinf"], nurange[0] - margin, side="right")
                - 1
            )
            # imax :
            # index i of the "w(i)-w(i+1)" element in dic_def["numtag"] such
            # that w(i+1)<=nurange[1]
            imax = (
                np.searchsorted(dic_def["numinf"], nurange[1] + margin, side="right")
                - 1
            )
            # not to exceed index out of the range
            imax = np.min([imax, len(dic_def["numinf"]) - 2])
            self.trans_file = []
            self.num_tag = []

            for k, i in enumerate(range(imin, imax + 1)):
                trans_file = self.path / pathlib.Path(
                    molec + "__" + dic_def["numtag"][i] + ".trans.bz2"
                )
                self.trans_file.append(trans_file)
                self.num_tag.append(dic_def["numtag"][i])

        # some verbose
        if self.verbose:
            print("Molecule: ", molecule)
            print("Isotopologue: ", self.isotope_fullname)
            print("ExoMol database: ", database)
            print("Local folder: ", self.path)
            print("Transition files: ")

        # Look-up missing parameters and write file
        # -----------------------------------------
        for trans_file, num_tag in zip(self.trans_file, self.num_tag):
            if self.verbose:
                print(
                    f"\t => File {os.path.splitext(os.path.basename((trans_file)))[0]}"
                )

            if cache == "regen" and mgr.cache_file(trans_file).exists():
                if self.verbose:
                    print("\t\t => `regen = True`. Removing the file")
                os.remove(mgr.cache_file(trans_file))

            if not mgr.cache_file(trans_file).exists():
                if cache == "force":
                    raise ValueError(
                        f"Cache file {str(mgr.cache_file(trans_file))} does not exist"
                    )

                if not trans_file.exists():
                    self.download(molec, extension=[".trans.bz2"], numtag=num_tag)
                if self.verbose:
                    print(
                        f"\t\t => Caching the *.trans.bz2 file to the {engine} (*.h5) format. After the second time, it will become much faster."
                    )
                    print(f"\t\t => You can deleted the 'trans.bz2' file by hand.")
                trans = read_trans(
                    trans_file, engine="vaex" if engine == "vaex" else "csv"
                )
                # TODO: add option to delete file at the end

                # Complete transition data with lookup on upper & lower state :
                # In particular, compute gup and elower

                trans = pickup_gE(
                    states,
                    trans,
                    dic_def,
                    skip_optional_data=skip_optional_data,
                    engine=engine,
                )

                ##Recompute Line strength:
                from radis.lbl.base import (  # TODO: move elsewhere
                    linestrength_from_Einstein,
                )

                self.Sij0 = linestrength_from_Einstein(
                    A=trans["A"],
                    gu=trans["gup"],
                    El=trans["elower"],
                    Ia=1,  #  Sij0 is a linestrength calculated without taking into account isotopic abundance (unlike line intensity parameter of HITRAN. In RADIS this is corrected for in fetch_exomol()  )
                    nu=trans["nu_lines"],
                    Q=self.QTref,
                    T=self.Tref,
                )

                trans["Sij0"] = self.Sij0

                mgr.write(mgr.cache_file(trans_file), trans)

    def set_broadening_coef(
        self,
        df,
        alpha_ref_def=None,
        n_Texp_def=None,
        output=None,
        add_columns=True,
        species=None,
    ):
        """setting broadening parameters

        Parameters
        ----------
        df: Data Frame
        alpha_ref: set default alpha_ref and apply it. None=use self.alpha_ref_def
        n_Texp_def: set default n_Texp and apply it. None=use self.n_Texp_def
        add_columns: adds alpha_ref and n_Texp columns to df
        species: to select which broadener will be used. Default is "air".

        Returns
        -------
        None. Store values in Data Frame.

        """
        if self.engine == "vaex":
            import vaex

        if alpha_ref_def:
            self.alpha_ref_def = alpha_ref_def
        if n_Texp_def:
            self.n_Texp_def = n_Texp_def

        file = self.broad_files[species]

        if self.broadf and os.path.exists(file):
            bdat = read_broad(file, output)

            if self.verbose > 1:
                print(f"The file `{os.path.basename(file)}` is used.")

            codelv = check_code_level(bdat, output=output)
            if self.verbose:
                print("Broadening code level:", codelv)

            if codelv == "a0":
                j2alpha_ref, j2n_Texp = make_j2b(
                    bdat,
                    alpha_ref_default=self.alpha_ref_def,
                    n_Texp_default=self.n_Texp_def,
                    jlower_max=df["jlower"].max(),
                    output=output,
                )
                self.alpha_ref = j2alpha_ref[df["jlower"].values.astype(int)]
                self.n_Texp = j2n_Texp[df["jlower"].values.astype(int)]
            elif codelv == "a1":
                j2alpha_ref, j2n_Texp = make_j2b(
                    bdat,
                    alpha_ref_default=self.alpha_ref_def,
                    n_Texp_default=self.n_Texp_def,
                    jlower_max=df["jlower"].max(),
                    output=output,
                )
                jj2alpha_ref, jj2n_Texp = make_jj2b(
                    bdat,
                    j2alpha_ref_def=j2alpha_ref,
                    j2n_Texp_def=j2n_Texp,
                    jupper_max=df["jupper"].max(),
                    output=output,
                )
                self.alpha_ref = (
                    np.array(
                        jj2alpha_ref[
                            df["jlower"].values.astype(int),
                            df["jupper"].values.astype(int),
                        ]
                    )
                    if self.engine != "vaex"
                    else jj2alpha_ref[
                        df["jlower"].values.astype(int), df["jupper"].values.astype(int)
                    ]
                )
                self.n_Texp = (
                    np.array(
                        jj2n_Texp[
                            df["jlower"].values.astype(int),
                            df["jupper"].values.astype(int),
                        ]
                    )
                    if self.engine != "vaex"
                    else jj2n_Texp[
                        df["jlower"].values.astype(int), df["jupper"].values.astype(int)
                    ]
                )
            elif codelv == "m0":

                # label P, Q, and R the transitions
                df["PQR"] = df["jupper"] - df["jlower"]  # P:+1, Q:0, R:-1

                if output == "vaex":
                    df["m"] = df.func.where(
                        df["PQR"] == -1,
                        -df["jlower"],
                        df.func.where(
                            df["PQR"] == 0,
                            df["jlower"],
                            df.func.where(df["PQR"] == 1, df["jlower"] + 1, 0),
                        ),
                    )
                    # np.nan causes the error in vaex when using map,
                    # 0 is regarded as the exeption value because m must not be zero
                else:
                    df["m"] = np.where(
                        df["PQR"] == -1,
                        -df["jlower"],
                        np.where(
                            df["PQR"] == 0,
                            df["jlower"],
                            np.where(df["PQR"] == 1, df["jlower"] + 1, np.nan),
                        ),
                    )

                    # Check for values outside -1, 0, 1 and raise a warning
                    invalid_pqr = df[~df["PQR"].isin([-1, 0, 1])]
                    # Not working for Vaex
                    if not invalid_pqr.empty:
                        warnings.warn(
                            f"Found {len(invalid_pqr)} values in 'PQR' outside of -1, 0, 1: {invalid_pqr['PQR'].unique()}"
                        )
                alpha_ref_dict = dict(zip(bdat["jlower"], bdat["alpha_ref"]))
                self.alpha_ref = (
                    np.array(df["m"].map(alpha_ref_dict).values)
                    if self.engine != "vaex"
                    else df["m"].map(alpha_ref_dict).values
                )

                n_Texp_dict = dict(zip(bdat["jlower"], bdat["n_Texp"]))
                self.n_Texp = (
                    np.array(df["m"].map(n_Texp_dict).values)
                    if self.engine != "vaex"
                    else df["m"].map(n_Texp_dict).values
                )
                ## for pandas but returns DataFrame
                # bdat.set_index("jlower", inplace=True)
                # self.alpha_ref = df["m"].map(bdat["alpha_ref"])
                # self.n_Texp = df["m"].map(bdat["n_Texp"])
                ## fill values outside of m range (but this gives DataFrame instead of np.array)
                # self.alpha_ref = self.alpha_ref.fillna(self.alpha_ref_def)
                # self.n_Texp = self.n_Texp.fillna(self.n_Texp_def)

            else:
                warnings.warn(
                    f"The broadening file contains this broadening code: {codelv}."
                    + " This broadening code is NOT implemented yet.\n"
                    + "Using default parameters instead."
                )
                if output != "vaex":
                    self.alpha_ref = self.alpha_ref_def * np.ones(len(df))
                    self.n_Texp = self.n_Texp_def * np.ones(len(df))
                else:
                    self.alpha_ref = vaex.from_arrays(
                        alpha_ref=[alpha_ref_def] * np.ones(len(df))
                    )
                    self.n_Texp = vaex.from_arrays(
                        n_Texp=[n_Texp_def] * np.ones(len(df))
                    )
        else:
            if not os.path.exists(file):
                warnings.warn(
                    f"Could not load `{os.path.basename(file)}`. The default broadening parameters are used.\n"
                )
            print("The default broadening parameters are used.")

            if output != "vaex":
                self.alpha_ref = np.array(self.alpha_ref_def * np.ones(len(df)))
                self.n_Texp = np.array(self.n_Texp_def * np.ones(len(df)))
            else:
                self.alpha_ref = vaex.from_arrays(
                    alpha_ref=[alpha_ref_def] * np.ones(len(df))
                )
                self.n_Texp = vaex.from_arrays(n_Texp=[n_Texp_def] * np.ones(len(df)))

        # Status: the 2 solumns self.alpha_ref and self.n_Texp are ready.
        # Next step: add the 2 solumns in df with the proper labels
        if add_columns:
            if species == "air":
                self.add_column(df, "airbrd", self.alpha_ref)
                self.add_column(df, "Tdpair", self.n_Texp)
            elif species == "self":
                self.add_column(df, "selbrd", self.alpha_ref)
                self.add_column(df, "selbrd_Tdpair", self.n_Texp)
            else:
                raise NotImplementedError(
                    "Please post on https://github.com/radis/radis to ask for this feature."
                )

    def QT_interp(self, T):
        """interpolated partition function

        Parameters
        ----------
        T: temperature

        Returns
        -------
        Q(T) interpolated in jnp.array

        """
        return np.interp(T, self.T_gQT, self.gQT)

    def qr_interp(self, T):
        """interpolated partition function ratio

        Parameters
        ----------
        T: temperature

        Returns
        -------
        qr(T)=Q(T)/Q(Tref) interpolated in jnp.array

        """
        return self.QT_interp(T) / self.QT_interp(self.Tref)

    def _calculate_download_size(self, url, pfname_arr):
        """Calculates the total download size of all files"""
        total_download_size_bytes = 0

        for pfname in pfname_arr:
            try:
                request = urllib.request.Request(url + pfname, method="HEAD")
                with urllib.request.urlopen(request) as response:
                    file_size = response.headers.get("Content-Length")

                    if file_size is None:
                        warnings.warn(
                            f"Missing Content-Length for {pfname}. Skipping.",
                            UserWarning,
                        )
                        continue

                    file_size = int(file_size)
                    total_download_size_bytes += file_size
            except Exception as e:
                warnings.warn(f"Failed to fetch size for {pfname}: {e}", UserWarning)

        total_download_size_gb = total_download_size_bytes / (1024**3)

        return total_download_size_gb

    def _construct_filenames_and_url(self, ext, molname_simple, tag, molec, numtag):
        if ext == ".trans.bz2" and numtag is not None:
            ext = "__" + numtag + ext

        if ext == ".broad":
            partners_success = np.ones(len(self.broad_partners), dtype=bool)
            pfname_arr = [tag[0] + "__" + s + ext for s in self.broad_partners]
            url = EXOMOL_URL + molname_simple + "/" + tag[0] + "/"
            return pfname_arr, url, partners_success
        else:
            pfname_arr = [molec + ext]
            url = EXOMOL_URL + molname_simple + "/" + tag[0] + "/" + tag[1] + "/"
            return pfname_arr, url

    def download(self, molec, extension, numtag=None):
        """Downloading Exomol files

        Parameters
        ----------
        molec: like "12C-16O__Li2015"
        extension: extension list e.g. [".pf",".def",".trans.bz2",".states.bz2",".broad"]
        numtag: number tag of transition file if exists. e.g. "11100-11200"

        Notes
        -----
        The download URL is written in exojax.utils.url.

        """
        import os
        import urllib.request

        tag = molec.split("__")
        molname_simple = exact_molname_exomol_to_simple_molname(tag[0])
        # TODO: add progress bar
        total_download_size_gb = 0
        for ext in extension:
            if ext == ".trans.bz2" and numtag is not None:
                ext = "__" + numtag + ext

            if ext == ".broad":
                pfname_arr, url, partners_success = self._construct_filenames_and_url(
                    ext, molname_simple, tag, molec, numtag
                )
            else:
                pfname_arr, url = self._construct_filenames_and_url(
                    ext, molname_simple, tag, molec, numtag
                )

            total_download_size_gb += self._calculate_download_size(url, pfname_arr)

        if self.verbose:
            print(
                f"Total download size {pfname_arr} is: {total_download_size_gb:.6f} GB"
            )

        from radis import config

        MAX_SIZE_GB = config["WARN_LARGE_DOWNLOAD_ABOVE_X_GB"]

        if total_download_size_gb > MAX_SIZE_GB:
            warning_msg = (
                f"The total download size is {total_download_size_gb:.2f} GB, which will take time and potential a significant portion of your disk memory."
                "To prevent this warning, you increase the limit using `radis.config['WARN_LARGE_DOWNLOAD_ABOVE_X_GB'] =  1`."
            )
            warnings.warn(warning_msg, UserWarning)

        for ext in extension:
            if ext == ".trans.bz2" and numtag is not None:
                ext = "__" + numtag + ext

            if ext == ".broad":
                pfname_arr, url, partners_success = self._construct_filenames_and_url(
                    ext, molname_simple, tag, molec, numtag
                )
            else:
                pfname_arr, url = self._construct_filenames_and_url(
                    ext, molname_simple, tag, molec, numtag
                )

            for index, pfname in enumerate(pfname_arr):
                pfpath = url + pfname
                os.makedirs(str(self.path), exist_ok=True)
                if self.verbose:
                    print(
                        f"\t\t => Downloading from {pfpath}"
                    )  # modify indent accordingly print in __init__
                try:
                    urllib.request.urlretrieve(pfpath, str(self.path / pfname))
                except HTTPError:
                    if ext == ".broad":
                        partners_success[index] = False
                    print(f"Error: Couldn't download {ext} file at {pfpath} and save.")

            if ext == ".broad" and self.verbose:
                patners_arr = np.array(self.broad_partners)
                print(
                    f"\nSummary of broadening files downloaded:\n\tSuccess: {patners_arr[partners_success]}\n\tFail: {patners_arr[~partners_success]}\n"
                )

    def to_partition_function_tabulator(self):
        """Generate a :py:class:`~radis.levels.partfunc.PartFuncExoMol` object"""
        from radis.levels.partfunc import PartFuncExoMol

        return PartFuncExoMol(self.isotope_fullname, self.T_gQT, self.gQT)


if __name__ == "__main__":
    print(exact_molname_exomol_to_simple_molname("12C-1H4"))
    print(exact_molname_exomol_to_simple_molname("23Na-16O-1H"))
    print(exact_molname_exomol_to_simple_molname("HeH_p"))
    print(exact_molname_exomol_to_simple_molname("trans-31P2-1H-2H"))  # not working

    # TODO : move to unitary tests of exomol_utils
    assert exact_molname_exomol_to_simple_molname("12C-1H4") == "CH4"
    assert exact_molname_exomol_to_simple_molname("23Na-16O-1H") == "NaOH"
    assert exact_molname_exomol_to_simple_molname("HeH_p") == "HeH_p"
    assert (
        exact_molname_exomol_to_simple_molname("trans-31P2-1H-2H") == "trans-31P2-1H-2H"
    )  # convert not working

    assert exact_molname_exomol_to_simple_molname("12C-16O2") == "CO2"
    assert exact_molname_exomol_to_simple_molname("13C-16O2") == "CO2"

    # mdb=MdbExomol("/home/kawahara/exojax/data/CO/12C-16O/Li2015/")
    # mdb=MdbExomol("/home/kawahara/exojax/data/CH4/12C-1H4/YT34to10/",nurange=[6050.0,6150.0])

    # %% Test  by overriding Spectrumfactory's DataFrame df0
    # mdb = MdbExomol(".database/H2O/1H2-16O/POKAZATEL", [4310.0, 4320.0], crit=1.0e-45)
    # df = mdb.to_df()
    # from radis import SpectrumFactory
    # sf = SpectrumFactory(
    #     4310,
    #     4320,
    #     molecule="H2O",
    #     isotope="1",
    # )
    # sf.fetch_databank(
    #     "hitran"
    # )  # placeholder. Load lines (will be replaced), load partition function.
    # sf.df0 = df  # override.
    # s = sf.eq_spectrum(500, name="ExoMol")
    # # sf.fetch_databank('hitran')  # placeholder. Load lines (will be replaced), load partition function.
    # # s_hit = sf.eq_spectrum(500, name='HITRAN')

    # %% Test by direct calculation
    # import pytest

    # print("Testing factory:", pytest.main(["../test/io/test_exomol.py"]))

    # %% RADIS-like Example
    # uses fetch_exomol() internally

    from radis import calc_spectrum

    s = calc_spectrum(
        wavelength_min=1.630e4,
        wavelength_max=1.6305e4,
        molecule="CH4",
        isotope="1",
        pressure=1.01325,  # bar
        Tgas=1000,  # K
        mole_fraction=0.1,
        path_length=1,  # cm
        databank=(
            "exomol",
            "YT10to10",
        ),  # Simply use 'exomol' for the recommended database
    )
    # s.apply_slit(1, "cm-1")  # simulate an experimental slit
    s.plot("xsection")

    # %% Exojax like Example

    class mdbExoMol:

        # hardcode attribute names, to prevent typos and the declaration of unwanted parameters
        __slots__ = [
            "Sij0",
            "logsij0",
            "nu_lines",
            "A",
            "elower",
            "eupper",
            "gupper",
            "jlower",
            "jupper",
        ]

        def __init__(
            self,
            molecule,
            path,
            nurange=[-np.inf, np.inf],
            local_databases="~/exojax",
        ):
            """
            Parameters
            ----------
            molecule: molecule name
            path : local path, mirror of ExoMol path
            nurange : TYPE, optional
                DESCRIPTION. The default is [-np.inf, np.inf].
            crit : TYPE, optional
                DESCRIPTION. The default is -np.inf.

            Returns
            -------
            DataFrame

            Examples
            --------
            ::

                mdbCH4 = mdbExoMol("CH4", '.database/CH4/12C-1H4/YT10to10/', nus, crit=1.e-30)
                print(len(mdbCH4.nu_lines), "lines")
                mdbCH4.elower

            Available columns::

                [
                    "Sij0",
                    "logsij0",
                    "nu_lines",
                    "A",
                    "elower",
                    "eupper",
                    "gupper",
                    "jlower",
                    "jupper",
                ]

            """

            wavenum_min, wavenum_max = np.min(nurange), np.max(nurange)
            if wavenum_min == -np.inf:
                wavenum_min = None
            if wavenum_max == np.inf:
                wavenum_max = None

            # Set-up database, download files and set-up cache files if needed
            mdb = MdbExomol(
                path,
                molecule=molecule,
                local_databases=local_databases,
                nurange=[wavenum_min, wavenum_max],
            )

            # Get cache files to load :
            mgr = mdb.get_datafile_manager()
            local_files = [mgr.cache_file(f) for f in mdb.trans_file]

            # Load them: #temperature dependent criterion
            jdict = mdb.load(
                local_files,
                columns=[k for k in self.__slots__ if k not in ["logsij0"]],
                lower_bound=([("nu_lines", wavenum_min)] if wavenum_min else [])
                + ([("Sij0", mdb.crit)] if not np.isneginf(mdb.crit) else []),
                upper_bound=([("nu_lines", wavenum_max)] if wavenum_max else []),
                output="jax",
            )

            # set attributes, accessible as e.g:  mdb.nu_lines
            for k in jdict.keys():
                setattr(self, k, jdict[k])

    nus = np.linspace(1e7 / 1.630e4, 1e7 / 1.6305e4)

    # Download new ExoMol repo (in ~/exomol)
    mdbCH4 = mdbExoMol(
        "CH4",
        ".database/CH4/12C-1H4/YT10to10/",
        nus,
        crit=1.0e-30,
        local_databases=".",  # use local folder
    )

    print(len(mdbCH4.nu_lines), "lines")
    mdbCH4.elower

    # Or use RADIS's folder  (# by default ~/.radisdb/exomol)
    import radis

    mdbCH4_2 = mdbExoMol(
        "CH4",
        "CH4/12C-1H4/YT10to10/",
        nus,
        crit=1.0e-30,
        local_databases=pathlib.Path(radis.config["DEFAULT_DOWNLOAD_PATH"]) / "exomol",
    )
    # ... ready to run Jax calculations

    # %%

    # """ExoMol lines can be downloaded and accessed separately using
    # :py:func:`~radis.io.exomol.fetch_exomol`
    # """

    # # See line data:
    # from radis.io.exomol import fetch_exomol

    # df = fetch_exomol("SiO", database="EBJT", isotope="1", load_wavenum_max=5000)
    # print(df)

    # #%%
    # # See the list of recommended databases for the 1st isotope of SiO :
    # from radis.api.exomolapi import get_exomol_database_list, get_exomol_full_isotope_name

    # databases, recommended = get_exomol_database_list(
    #     "SiO", get_exomol_full_isotope_name("SiO", 1)
    # )
    # print("Databases for SiO: ", databases)
    # print("Database recommended by ExoMol: ", recommended)
