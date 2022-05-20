"""API for Exomol molecular database

Borrowed from the `Exojax <https://github.com/HajimeKawahara/exojax>`__
code (which you should also have a look at !), by @HajimeKawahara, under MIT License.

"""
import bz2
import re

import numpy as np
import pandas as pd

try:
    from .hdf5 import vaexsafe_colname
except:
    from radis.io.hdf5 import vaexsafe_colname

from radis.misc.warning import InconsistentDatabaseError


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

    # SOME DEF FILES CONTAINS ERRORS. THESE ARE THE EXCEPTIONS
    if deff.stem == "12C-16O2__UCL-4000":
        ntransf = 20
    if deff.stem == "14N-1H3__CoYuTe":
        maxnu = 20000.0
    if deff.stem == "12C2-1H2__aCeTY":
        if molmass == 12.0:
            molmass = 26.0
            print(
                f"Known error in ExoMol def file, molmass corrected from 12.0 to {molmass}"
            )
        if quantum_labels == [
            "totalSym",
            "v1",
            "v2",
            "v3",
            "v4",
            "v5",
            "v5",
            "v7",
            "vibSym",
            "K",
            "rotSym",
        ]:
            quantum_labels = [
                "totalSym",
                "v1",
                "v2",
                "v3",
                "v4",
                "v5",
                "v6",
                "v7",
                "vibSym",
                "K",
                "rotSym",
            ]
            print(
                f"Known error in ExoMol def file, quantum_labels corrected from '['totalSym', 'v1', 'v2', 'v3', 'v4', 'v5', 'v5', 'v7', 'vibSym', 'K', 'rotSym']' to {quantum_labels}"
            )

    if ntransf > 1:
        dnufile = maxnu / ntransf
        numinf = dnufile * np.array(range(ntransf + 1))
        numtag = []
        for i in range(len(numinf) - 1):
            imin = "{:05}".format(int(numinf[i]))
            imax = "{:05}".format(int(numinf[i + 1]))
            numtag.append(imin + "-" + imax)
    else:
        numinf = None
        numtag = ""

    output = {
        "n_Texp": n_Texp,
        "alpha_ref": alpha_ref,
        "molmass": molmass,
        "numinf": numinf,
        "numtag": numtag,
        "quantum_labels": quantum_labels,  # array
        "Landé": lande,  # bool
        "lifetime": lifetime,  # bool
    }
    return output


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
                    f"Error reading {transf}. Maybe the file was corrupted during download ? You can try to delete it."
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
                    f"Error reading {transf}. Maybe the file was corrupted during download ? You can try to delete it."
                ) from err
    else:
        raise NotImplementedError(engine)

    return dat


def read_states(statesf, dic_def, engine="vaex", skip_optional_data=True):
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
    N = np.size(dic_def["quantum_labels"])
    quantum_labels = dic_def["quantum_labels"]

    mandatory_usecol = np.arange(4)
    mandatory_fields = ("i", "E", "g", "J")
    if skip_optional_data:
        usecol = mandatory_usecol
        names = mandatory_fields
    else:
        if dic_def["Landé"] and dic_def["lifetime"]:
            usecol = np.concatenate((mandatory_usecol, 5 + 2 + np.arange(N)))
            names = mandatory_fields + ("Lande", "lifetime") + tuple(quantum_labels)
        elif dic_def["Landé"]:
            usecol = np.concatenate((mandatory_usecol, 5 + 1 + np.arange(N)))
            names = mandatory_fields + ("Lande",) + tuple(quantum_labels)
        elif dic_def["lifetime"]:
            usecol = np.concatenate((mandatory_usecol, 5 + 1 + np.arange(N)))
            names = mandatory_fields + ("lifetime",) + tuple(quantum_labels)
        else:  # no lifetime, nor landé according to the def file
            usecol = np.concatenate((mandatory_usecol, 5 + np.arange(N)))
            names = mandatory_fields + tuple(quantum_labels)
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
            raise InconsistentDatabaseError(
                "The EXOMOL definitions and states files are inconsistent.\n"
                + "Some data are specified as available by the *.def file, but are absent in the *.bz2 file.\n"
                + "Set `skip_optional_data=False` in `fetch_exomol()` to load only the required data for a LTE computation.\n"
                + "The problematic optional data/columns will be ignored."
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
                    f"Error reading {statesf}. Maybe the file was corrupted during download ? You can try to delete it."
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

    return dat


def pickup_gE(
    states, trans, trans_file, dic_def, trans_lines=False, skip_optional_data=True
):
    """extract g_upper (gup), E_lower (elower), and J_lower and J_upper from states
    DataFrame and insert them into the transition DataFrame.

    Parameters
    ----------
    states: states DataFrame  - the i, E, g, J are in the 4 first columns
    trans: transition numpy array
    trans_file: name of the transition file
    trans_lines: By default (False) we use nu_lines computed using the state file, i.e. E_upper - E_lower. If trans_nuline=True, we use the nu_lines in the transition file. Note that some trans files do not this info.
    dic_def: Informations about additional quantum labels
    add_quantum_labels: bool . If True fetch all quantum labels in dic_def['quantum_labels'] from states into transitions

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

    def map_add(col, new_col, trans_key, states_key="i"):
        """Lookup `key` in states and add it in trans, using the level ``i``
        for upper and lower state

        Examples
        --------
        ::

            map_add("E", "E_lower", "i_lower")
        """
        try:  # pytable
            trans[new_col] = trans[trans_key].map(dict(states[col]))
        except:  # a priori, vaex version  (TODO : replace with dict() approach in vaex too)
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

    map_add("g", "gup", "i_upper")
    map_add("J", "jlower", "i_lower")
    map_add("J", "jupper", "i_upper")

    map_add("E", "elower", "i_lower")

    def has_nan(column):
        try:  # Vaex
            return column.countnan() > 0
        except AttributeError:  # Pandas
            return column.hasnans

    if not "nu_lines" in trans or has_nan(trans["nu_lines"]):
        map_add("E", "eupper", "i_upper")
        trans["nu_lines"] = trans["eupper"] - trans["elower"]

    ### Step 2. Extra quantum numbers (e/f parity, vib and rot numbers)
    # -----------------------------------------------------------------
    if not skip_optional_data:
        for q in dic_def["quantum_labels"]:
            map_add(q, f"{q}_l", "i_lower")
            # map_add(q, f"{q}_u", "i_upper")

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


def read_broad(broadf):
    """Reading braodening file (.broad)
    Parameters
    ----------
    broadf: .broad file

    Returns
    -------
    broadening info in bdat form (pandas), defined by this instance.

    Notes
    -----
    See Table 16 in https://arxiv.org/pdf/1603.05890.pdf
    """
    bdat = pd.read_csv(
        broadf,
        sep=r"\s+",
        names=(
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
        ),
    )

    return bdat


def check_bdat(bdat):
    """cheking codes in .broad
    Args:
       bdat: exomol .broad data given by exomolapi.read_broad
    Returns:
       code level: None, a0, a1, other codes unavailable currently,
    """

    def checkcode(code):
        cmask = bdat["code"] == code
        if len(bdat["code"][cmask]) > 0:
            return True
        else:
            return False

    codelv = None
    for code in ["a0", "a1"]:
        if checkcode(code):
            codelv = code

    return codelv


def make_j2b(bdat, alpha_ref_default=0.07, n_Texp_default=0.5, jlower_max=None):
    """compute j2b (code a0, map from jlower to alpha_ref)

    Args:
       bdat: exomol .broad data given by exomolapi.read_broad
       alpha_ref_default: default value
       n_Texp_default: default value
       jlower_max: maximum number of jlower
    Returns:
       j2alpha_ref[jlower] provides alpha_ref for jlower
       j2n_Texp[jlower]  provides nT_exp for jlower
    """
    # a0
    cmask = bdat["code"] == "a0"
    jlower_arr = np.array(bdat["jlower"][cmask], dtype=int)
    alpha_ref_arr = np.array(bdat["alpha_ref"][cmask])
    n_Texp_arr = np.array(bdat["n_Texp"][cmask])

    if jlower_max is None:
        Nblower = np.max(jlower_arr) + 1
    else:
        Nblower = np.max([jlower_max, np.max(jlower_arr)]) + 1
    j2alpha_ref = np.ones(Nblower) * alpha_ref_default
    j2n_Texp = np.ones(Nblower) * n_Texp_default

    j2alpha_ref[jlower_arr] = alpha_ref_arr
    j2n_Texp[jlower_arr] = n_Texp_arr

    Ndef = Nblower - (np.max(jlower_arr) + 1)
    if Ndef > 0:
        print(
            "default broadening parameters are used for ",
            Ndef,
            " J lower states in ",
            Nblower,
            " states",
        )

    return j2alpha_ref, j2n_Texp


def make_jj2b(bdat, j2alpha_ref_def, j2n_Texp_def, jupper_max=None):
    """compute jj2b (code a1, map from (jlower, jupper) to alpha_ref and n_Texp)

    Args:
       bdat: exomol .broad data given by exomolapi.read_broad
       j2alpha_ref_def: default value from a0
       j2n_Texp_def: default value from a0
       jupper_max: maximum number of jupper
    Returns:
       jj2alpha_ref[jlower,jupper] provides alpha_ref for (jlower, jupper)
       jj2n_Texp[jlower,jupper]  provides nT_exp for (jlower, jupper)
    Note:
       The pair of (jlower, jupper) for which broadening parameters are not given, jj2XXX contains None.
    """
    # a1
    cmask = bdat["code"] == "a1"
    jlower_arr = np.array(bdat["jlower"][cmask], dtype=int)
    jupper_arr = np.array(bdat["jupper"][cmask], dtype=int)
    alpha_ref_arr = np.array(bdat["alpha_ref"][cmask])
    n_Texp_arr = np.array(bdat["n_Texp"][cmask])

    if jupper_max is None:
        Nbupper = np.max(jupper_arr) + 1
    else:
        Nbupper = np.max([jupper_max, np.max(jupper_arr)]) + 1

    jj2alpha_ref = j2alpha_ref_def[:, np.newaxis] * np.ones(Nbupper)
    jj2n_Texp = j2n_Texp_def[:, np.newaxis] * np.ones(Nbupper)

    jj2alpha_ref[jlower_arr, jupper_arr] = alpha_ref_arr
    jj2n_Texp[jlower_arr, jupper_arr] = n_Texp_arr

    return jj2alpha_ref, jj2n_Texp


if __name__ == "__main__":
    pass
