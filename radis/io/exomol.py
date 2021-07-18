"""Molecular database (MDB) class

   * MdbExomol is the MDB for ExoMol

Initial code borrowed from the `Exojax <https://github.com/HajimeKawahara/exojax>`__
code (which you should also have a look at !), by @HajimeKawahara, under MIT License.

"""
import pathlib
from ntpath import join
from os.path import expanduser

import numpy as np

# from exojax.spec import hapi, exomolapi, exomol
# from exojax.spec.hitran import gamma_natural as gn
import pandas as pd

try:
    from . import exomolapi
    from .exomol_utils import e2s
except:  # if local import
    from radis.io import exomolapi
    from radis.io.exomol_utils import e2s

from radis.db.classes import get_molecule_identifier

EXOMOL_URL = "http://www.exomol.com/db/"

KNOWN_EXOMOL_ISOTOPES_NAMES = {
    "H2O": {1: "1H2-16O"},
    "NH3": {1: "14N-1H3"},
    "H2S": {1: "1H2-32S"},
    "FeH": {1: "56Fe-1H"},
}
"""molecule: isotope: full name"""
KNOWN_EXOMOL_DATABASE_NAMES = {
    "H2O": {"1H2-16O": ["POKAZATEL"]},
    "NH3": {"14N-1H3": ["CoYuTe"]},
    "H2S": {"1H2-32S": ["AYT2"]},
    "FeH": {"56Fe-1H": ["MoLLIST"]},
}


def fetch_exomol(
    molecule,
    database=None,
    local_databases="~/.radisdb/exomol/",
    databank_name="EXOMOL-{molecule}",
    isotope=None,
    load_wavenum_min=None,
    load_wavenum_max=None,
    cache=True,
    verbose=True,
    clean_cache_files=True,
    return_local_path=False,
):
    """Stream ExoMol file from EXOMOL website. Unzip and build a HDF5 file directly.

    Returns a Pandas DataFrame containing all lines.

    Parameters
    ----------
    molecule: str
        ExoMol molecule
    database: str
        database name. Ex:: ``POKAZATEL`` or ``BT2`` for ``H2O``. See
        :py:data:`~radis.io.exomol.KNOWN_EXOMOL_DATABASE_NAMES`. If ``None`` and
        there is only one database available, use it.
    local_databases: str
        where to create the RADIS HDF5 files. Default ``"~/.radisdb/exomol"``
    databank_name: str
        name of the databank in RADIS :ref:`Configuration file <label_lbl_config_file>`
        Default ``"EXOMOL-{molecule}"``
    isotope: str
        load only certain isotopes : ``'2'``, ``'1,2'``, etc. Default ``1``.
    load_wavenum_min, load_wavenum_max: float (cm-1)
        load only specific wavenumbers.

    Other Parameters
    ----------------
    cache: bool, or ``'regen'``
        if ``True``, use existing HDF5 file. If ``False`` or ``'regen'``, rebuild it.
        Default ``True``.
    verbose: bool
    clean_cache_files: bool
        if ``True`` clean downloaded cache files after HDF5 are created.
    return_local_path: bool
        if ``True``, also returns the path of the local database file.

    Returns
    -------
    df: pd.DataFrame
        Line list
        A HDF5 file is also created in ``local_databases`` and referenced
        in the :ref:`RADIS config file <label_lbl_config_file>` with name
        ``databank_name``
    local_path: str
        path of local database file if ``return_local_path``

    Notes
    -----
    if using ``load_only_wavenum_above/below`` or ``isotope``, the whole
    database is anyway downloaded and uncompressed to ``local_databases``
    fast access .HDF5 files (which will take a long time on first call). Only
    the expected wavenumber range & isotopes are returned. The .HFD5 parsing uses
    :py:func:`~radis.io.hdf5.hdf2df`

    See Also
    --------
    :py:func:`~radis.io.hdf5.hdf2df`

    """
    full_molecule_name = KNOWN_EXOMOL_ISOTOPES_NAMES[molecule][isotope]
    known_exomol_databases = KNOWN_EXOMOL_DATABASE_NAMES[molecule][full_molecule_name]

    if database is None:
        if len(known_exomol_databases) == 1:
            database = known_exomol_databases[0]
            if verbose:
                print(f"Using ExoMol database {database} for {full_molecule_name}")
        else:
            raise KeyError(
                f"Choose one of the several databases available for {full_molecule_name} in ExoMol: {known_exomol_databases}"
            )
    else:
        if database not in KNOWN_EXOMOL_DATABASE_NAMES[molecule][full_molecule_name]:
            raise KeyError(
                f"{database} is not of the known available ExoMol databases for {full_molecule_name}. Choose one of : {known_exomol_databases}"
            )

    local_path = join(
        expanduser(local_databases), molecule, full_molecule_name, database
    )

    mdb = MdbExomol(local_path, [load_wavenum_min, load_wavenum_max])

    df = mdb.to_df()

    return (df, mdb.path) if return_local_path else df


class MdbExomol(object):
    """molecular database of ExoMol

    MdbExomol is a class for ExoMol.

    Attributes:
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

    """

    def __init__(
        self,
        path,
        nurange=[-np.inf, np.inf],
        margin=1.0,
        crit=-np.inf,
        bkgdatm="H2",
        broadf=True,
    ):
        """Molecular database for Exomol form

        Args:
           path: path for Exomol data directory/tag. For instance, "/home/CO/12C-16O/Li2015"
           nurange: wavenumber range list (cm-1) or wavenumber array
           margin: margin for nurange (cm-1)
           crit: line strength lower limit for extraction
           bkgdatm: background atmosphere for broadening. e.g. H2, He,
           broadf: if False, the default broadening parameters in .def file is used

        Note:
           The trans/states files can be very large. For the first time to read it, we convert it to the feather-format. After the second-time, we use the feather format instead.

        """
        explanation = "Note: Couldn't find the feather format. We convert data to the feather format. After the second time, it will become much faster."

        self.path = pathlib.Path(path)
        t0 = self.path.parents[0].stem
        molec = t0 + "__" + str(self.path.stem)
        self.bkgdatm = bkgdatm
        print("Background atmosphere: ", self.bkgdatm)
        molecbroad = t0 + "__" + self.bkgdatm

        self.crit = crit
        self.margin = margin
        self.nurange = [np.min(nurange), np.max(nurange)]
        self.broadf = broadf
        # Where exomol files are
        self.states_file = self.path / pathlib.Path(molec + ".states.bz2")
        self.pf_file = self.path / pathlib.Path(molec + ".pf")
        self.def_file = self.path / pathlib.Path(molec + ".def")
        self.broad_file = self.path / pathlib.Path(molecbroad + ".broad")

        if not self.def_file.exists():
            self.download(molec, extension=[".def"])
        if not self.pf_file.exists():
            self.download(molec, extension=[".pf"])
        if not self.states_file.exists():
            self.download(molec, extension=[".states.bz2"])
        if not self.broad_file.exists():
            self.download(molec, extension=[".broad"])

        # Add molecule name
        tag = molec.split("__")
        self.molecule = e2s(tag[0])
        self.isotope = 1  # Placeholder. TODO : impement parsing of other isotopes.

        # load def
        (
            self.n_Texp_def,
            self.alpha_ref_def,
            self.molmass,
            numinf,
            numtag,
        ) = exomolapi.read_def(self.def_file)
        #  default n_Texp value if not given
        if self.n_Texp_def is None:
            self.n_Texp_def = 0.5
        #  default alpha_ref value if not given
        if self.alpha_ref_def is None:
            self.alpha_ref_def = 0.07

        # load states
        if self.states_file.with_suffix(".feather").exists():
            states = pd.read_feather(self.states_file.with_suffix(".feather"))
        else:
            print(explanation)
            states = exomolapi.read_states(self.states_file)
            states.to_feather(self.states_file.with_suffix(".feather"))
        # load pf
        pf = exomolapi.read_pf(self.pf_file)
        self.gQT = pf["QT"].to_numpy()  # grid QT
        self.T_gQT = pf["T"].to_numpy()  # T forgrid QT

        # trans file(s)
        print("Reading transition file")
        if numinf is None:
            self.trans_file = self.path / pathlib.Path(molec + ".trans.bz2")
            if not self.trans_file.exists():
                self.download(molec, [".trans.bz2"])

            if self.trans_file.with_suffix(".feather").exists():
                trans = pd.read_feather(self.trans_file.with_suffix(".feather"))
            else:
                print(explanation)
                trans = exomolapi.read_trans(self.trans_file)
                trans.to_feather(self.trans_file.with_suffix(".feather"))
            # compute gup and elower
            (
                self._A,
                self.nu_lines,
                self._elower,
                self._gpp,
                self._jlower,
                self._jupper,
            ) = exomolapi.pickup_gE(states, trans)
        else:
            imin = np.searchsorted(numinf, nurange[0], side="right") - 1  # left side
            imax = np.searchsorted(numinf, nurange[1], side="right") - 1  # left side
            self.trans_file = []
            for k, i in enumerate(range(imin, imax + 1)):
                trans_file = self.path / pathlib.Path(
                    molec + "__" + numtag[i] + ".trans.bz2"
                )
                if not trans_file.exists():
                    self.download(molec, extension=[".trans.bz2"], numtag=numtag[i])
                if trans_file.with_suffix(".feather").exists():
                    trans = pd.read_feather(trans_file.with_suffix(".feather"))
                else:
                    print(explanation)
                    trans = exomolapi.read_trans(trans_file)
                    trans.to_feather(trans_file.with_suffix(".feather"))
                self.trans_file.append(trans_file)
                # compute gup and elower
                if k == 0:
                    (
                        self._A,
                        self.nu_lines,
                        self._elower,
                        self._gpp,
                        self._jlower,
                        self._jupper,
                    ) = exomolapi.pickup_gE(states, trans)
                else:
                    Ax, nulx, elowerx, gppx, jlowerx, jupperx = exomolapi.pickup_gE(
                        states, trans
                    )
                    self._A = np.hstack([self._A, Ax])
                    self.nu_lines = np.hstack([self.nu_lines, nulx])
                    self._elower = np.hstack([self._elower, elowerx])
                    self._gpp = np.hstack([self._gpp, gppx])
                    self._jlower = np.hstack([self._jlower, jlowerx])
                    self._jupper = np.hstack([self._jupper, jupperx])

        self.Tref = 296.0
        self.QTref = np.array(self.QT_interp(self.Tref))

        Ia = 1  #  TODO    Add isotope abundance

        from radis.lbl.base import linestrength_from_Einstein  # TODO: move elsewhere

        self.Sij0 = linestrength_from_Einstein(
            A=self._A,
            gu=self._gpp,
            El=self._elower,
            Ia=Ia,
            nu=self.nu_lines,
            Q=self.QTref,
            T=self.Tref,
        )

        ### MASKING ###
        mask = (
            (self.nu_lines > self.nurange[0] - self.margin)
            * (self.nu_lines < self.nurange[1] + self.margin)
            * (self.Sij0 > self.crit)
        )

        self.masking(mask)

    def masking(self, mask):
        """applying mask and (re)generate jnp.arrays

        Args:
           mask: mask to be applied. self.mask is updated.

        """
        # TODO : replace with HDF5 masking ?

        # numpy float 64 Do not convert them jnp array
        self.nu_lines = self.nu_lines[mask]
        self.Sij0 = self.Sij0[mask]
        self._A = self._A[mask]
        self._elower = self._elower[mask]
        self._gpp = self._gpp[mask]
        self._jlower = self._jlower[mask]
        self._jupper = self._jupper[mask]

        # jnp arraysgamma_
        self.dev_nu_lines = np.array(self.nu_lines)
        self.logsij0 = np.array(np.log(self.Sij0))
        self.A = np.array(self._A)
        # self.gamma_natural = gn(self.A)   #  Natural Broadening neglected in RADIS
        self.elower = np.array(self._elower)
        self.gpp = np.array(self._gpp)
        self.jlower = np.array(self._jlower, dtype=int)
        self.jupper = np.array(self._jupper, dtype=int)
        ##Broadening parameters
        self.set_broadening()

    def set_broadening(self, alpha_ref_def=None, n_Texp_def=None):
        """setting broadening parameters

        Args:
           alpha_ref: set default alpha_ref and apply it. None=use self.alpha_ref_def
           n_Texp_def: set default n_Texp and apply it. None=use self.n_Texp_def
        """
        if alpha_ref_def:
            self.alpha_ref_def = alpha_ref_def
        if n_Texp_def:
            self.n_Texp_def = n_Texp_def

        if self.broadf:
            try:
                print(".broad is used.")
                bdat = exomolapi.read_broad(self.broad_file)
                codelv = exomolapi.check_bdat(bdat)
                print("Broadening code level=", codelv)
                if codelv == "a0":
                    j2alpha_ref, j2n_Texp = exomolapi.make_j2b(
                        bdat,
                        alpha_ref_default=self.alpha_ref_def,
                        n_Texp_default=self.n_Texp_def,
                        jlower_max=np.max(self._jlower),
                    )
                    self.alpha_ref = np.array(j2alpha_ref[self._jlower])
                    self.n_Texp = np.array(j2n_Texp[self._jlower])
                elif codelv == "a1":
                    j2alpha_ref, j2n_Texp = exomolapi.make_j2b(
                        bdat,
                        alpha_ref_default=self.alpha_ref_def,
                        n_Texp_default=self.n_Texp_def,
                        jlower_max=np.max(self._jlower),
                    )
                    jj2alpha_ref, jj2n_Texp = exomolapi.make_jj2b(
                        bdat,
                        j2alpha_ref_def=j2alpha_ref,
                        j2n_Texp_def=j2n_Texp,
                        jupper_max=np.max(self._jupper),
                    )
                    self.alpha_ref = np.array(jj2alpha_ref[self._jlower, self._jupper])
                    self.n_Texp = np.array(jj2n_Texp[self._jlower, self._jupper])
            except:
                print(
                    "Warning: Cannot load .broad. The default broadening parameters are used."
                )
                self.alpha_ref = np.array(
                    self.alpha_ref_def * np.ones_like(self._jlower)
                )
                self.n_Texp = np.array(self.n_Texp_def * np.ones_like(self._jlower))

        else:
            print("The default broadening parameters are used.")
            self.alpha_ref = np.array(self.alpha_ref_def * np.ones_like(self._jlower))
            self.n_Texp = np.array(self.n_Texp_def * np.ones_like(self._jlower))

    def QT_interp(self, T):
        """interpolated partition function

        Args:
           T: temperature

        Returns:
           Q(T) interpolated in jnp.array

        """
        return np.interp(T, self.T_gQT, self.gQT)

    def qr_interp(self, T):
        """interpolated partition function ratio

        Args:
           T: temperature

        Returns:
           qr(T)=Q(T)/Q(Tref) interpolated in jnp.array

        """
        return self.QT_interp(T) / self.QT_interp(self.Tref)

    def download(self, molec, extension, numtag=None):
        """Downloading Exomol files

        Args:
           molec: like "12C-16O__Li2015"
           extension: extension list e.g. [".pf",".def",".trans.bz2",".states.bz2",".broad"]
           numtag: number tag of transition file if exists. e.g. "11100-11200"

        Note:
           The download URL is written in exojax.utils.url.

        """
        import os
        import urllib.request

        tag = molec.split("__")
        molname_simple = e2s(tag[0])

        for ext in extension:
            if ext == ".trans.bz2" and numtag is not None:
                ext = "__" + numtag + ext

            if ext == ".broad":
                pfname_arr = [
                    tag[0] + "__H2" + ext,
                    tag[0] + "__He" + ext,
                    tag[0] + "__air" + ext,
                ]
                url = EXOMOL_URL + molname_simple + "/" + tag[0] + "/"
            else:
                pfname_arr = [molec + ext]
                url = EXOMOL_URL + molname_simple + "/" + tag[0] + "/" + tag[1] + "/"

            for pfname in pfname_arr:
                pfpath = url + pfname
                os.makedirs(str(self.path), exist_ok=True)
                print("Downloading " + pfpath)
                try:
                    urllib.request.urlretrieve(pfpath, str(self.path / pfname))
                except:
                    print("Error: Couldn't download " + ext + " file and save.")

    def to_df(self):
        """Export Line Database to a RADIS-friendly Pandas DataFrame"""

        ##Broadening parameters
        self.set_broadening()

        # TODO : define RADIS database format somewhere else; with description of the column names.
        df = pd.DataFrame(
            {
                "wav": self.nu_lines,
                "int": self.Sij0,
                "A": self._A,
                "airbrd": self.alpha_ref,  # temperature dependance exponent for air broadening
                # "selbrd": None,                  # self-temperature dependant exponent. Not implementedi in ExoMol
                "El": self._elower,
                # "Tdpair": None,
                # "Pshft": None,
                "ju": self._jupper,
                "gpp": self._gpp,
                "Tdpair": self.n_Texp,  # temperature dependance exponent. Here we use Tdair; no Tdpsel. TODO
            }
        )

        df.attrs["id"] = get_molecule_identifier(self.molecule)
        df.attrs["iso"] = 1  # TODO : get isotope number while parsing ExoMol name

        assert self.Tref == 296  # default in RADIS

        return df


if __name__ == "__main__":
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

    #%% Test by direct caclulation
    from radis import SpectrumFactory

    sf = SpectrumFactory(
        4310,
        4320,
        molecule="H2O",
        isotope="1",
        verbose=3,
    )
    sf.fetch_databank("exomol")
    s = sf.eq_spectrum(500, name="ExoMol")
    s.plot()

#    mask=mdb.A>1.e-42
#    mdb.masking(mask)
#    mdb=MdbExomol("/home/kawahara/exojax/data/exomol/NH3/14N-1H3/CoYuTe/",nurange=[6050.0,6150.0])
#    mdb=MdbExomol("/home/kawahara/exojax/data/exomol/H2S/1H2-32S/AYT2/",nurange=[6050.0,6150.0])
#    mdb=MdbExomol("/home/kawahara/exojax/data/exomol/FeH/56Fe-1H/MoLLIST/",nurange=[6050.0,6150.0])
#    mdb=MdbExomol("/home/kawahara/exojax/data/exomol/NO/14N-16O/NOname/14N-16O__NOname")
