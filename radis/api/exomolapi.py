"""Molecular database (MDB) class with multi-species broadening support"""

import os
import pathlib
import re
import warnings
import numpy as np
from urllib.request import urlopen, HTTPError
from bs4 import BeautifulSoup
from radis.api.dbmanager import DatabaseManager, get_auto_MEMORY_MAPPING_ENGINE
from radis.db.classes import EXOMOL_MOLECULES, EXOMOL_ONLY_ISOTOPES_NAMES

EXOMOL_URL = "https://www.exomol.com/db/"

# Utility Functions
def exact_molname_exomol_to_simple_molname(exact_exomol_molecule_name):
    """Convert ExoMol exact name (e.g., '12C-1H4') to simple formula ('CH4')"""
    try:
        molname_simple = ""
        for ele in exact_exomol_molecule_name.split("-"):
            alp = "".join(re.findall(r"\D", ele))
            num = re.split("[A-Z]", ele)[1] if re.split("[A-Z]", ele)[1].isdigit() else ""
            molname_simple += alp + num
        return "H2O" if molname_simple == "HHO" else molname_simple
    except Exception:
        warnings.warn(f"Could not convert {exact_exomol_molecule_name}")
        return exact_exomol_molecule_name

def get_exomol_full_isotope_name(molecule, isotope):
    """Get ExoMol isotope name from molecule and isotope number"""
    if (molecule, isotope) in EXOMOL_ONLY_ISOTOPES_NAMES:
        return EXOMOL_ONLY_ISOTOPES_NAMES[(molecule, isotope)]
    raise ValueError(f"Isotope not found for {molecule}({isotope})")

def get_list_of_known_isotopes(molecule):
    """List all isotopes for a molecule"""
    isotopes = []
    i = 1
    while True:
        try:
            isotopes.append(get_exomol_full_isotope_name(molecule, i))
            i += 1
        except ValueError:
            break
    return isotopes

def get_exomol_database_list(molecule, isotope_full_name=None):
    """Fetch available ExoMol databases for Molecule-Isotope pair"""
    if isotope_full_name is None:
        raise ValueError(f"Isotope name required. Known: {get_list_of_known_isotopes(molecule)}")
    
    try:
        url = f"{EXOMOL_URL}{molecule}/{isotope_full_name}"
        soup = BeautifulSoup(urlopen(url).read(), "lxml")
        databases = [a["title"] for a in soup.find_all("a", class_="list-group-item")]
        recommended = [a["title"] for a in soup.find_all("a", class_="recommended")]
        return list(set(databases)), recommended[0] if recommended else False
    except HTTPError as e:
        raise ValueError(f"Failed to fetch databases: {e}")

def parse_broad_file(broad_file_path):
    """Parse ExoMol .broad file into broadening coefficients"""
    broad_data = {}
    with open(broad_file_path, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                fields = line.split()
                if len(fields) >= 3:  # Assuming format: J, gamma, n
                    J, gamma, n = float(fields[0]), float(fields[1]), float(fields[2])
                    broad_data[J] = {'gamma': gamma, 'n': n}
    return broad_data

class MdbExomol(DatabaseManager):
    """Molecular database of ExoMol with multi-species broadening support"""
    
    DEFAULT_BROADENING = {'gamma': 0.1, 'n': 0.5}  # Default values if missing

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
        bkgdatm="Air",
        broadf=True,
        broadf_download=True,
        engine="vaex",
        verbose=True,
        cache=True,
        skip_optional_data=True,
        broadening_partners=None,  # Dictionary of diluents and mole fractions
    ):
        super().__init__(name, molecule, local_databases, engine, verbose=verbose)
        assert cache
        
        # Initialize paths
        self.path = pathlib.Path(path)
        if local_databases is not None:
            self.path = pathlib.Path(local_databases).expanduser() / self.path

        # Setup broadening partners
        self.broadening_partners = self._validate_broadening_partners(
            broadening_partners if broadening_partners is not None else {"air": 1.0}
        )
        self.broadening_coefficients = {}
        self.broadf = broadf
        self.broadf_download = broadf_download

        # Rest of initialization
        t0 = self.path.parents[0].stem
        molec = t0 + "__" + str(self.path.stem)
        self.crit = crit
        self.margin = margin
        self.nurange = [np.min(nurange), np.max(nurange)]
        self.wmin, self.wmax = np.min(nurange), np.max(nurange)
        self.engine = engine

        # Initialize file paths
        self._init_file_paths(molec, t0)
        
        # Load definition and partition function files
        self._load_def_and_pf_files(molec)
        
        # Load states and transition files
        self._load_states_and_transitions(molec, skip_optional_data)
        
        # Load broadening data
        if self.broadf:
            self._load_broadening_data(molec)

    def _validate_broadening_partners(self, broadening_partners):
        """Validate broadening partners dictionary and ensure mole fractions sum to 1"""
        if not isinstance(broadening_partners, dict):
            raise ValueError("broadening_partners must be a dict, e.g., {'air': 0.7, 'H2': 0.3}")
        total_fraction = sum(broadening_partners.values())
        if not np.isclose(total_fraction, 1.0, atol=1e-2):
            warnings.warn(f"Mole fractions sum to {total_fraction}, not 1. Normalizing.")
            for diluent in broadening_partners:
                broadening_partners[diluent] /= total_fraction
        return broadening_partners

    def _load_broadening_data(self, molec):
        """Load broadening coefficients for all specified diluents"""
        local_path = self.path.parents[0]
        for diluent in self.broadening_partners:
            broad_file = local_path / f"{molec}__{diluent}.broad"
            if broad_file.exists():
                self.broadening_coefficients[diluent] = parse_broad_file(broad_file)
            elif self.broadf_download:
                try:
                    self._download_broadening_file(molec, diluent, local_path)
                    self.broadening_coefficients[diluent] = parse_broad_file(broad_file)
                except Exception as e:
                    warnings.warn(f"Failed to download {diluent} broadening for {molec}: {e}. Using default.")
                    self.broadening_coefficients[diluent] = self.DEFAULT_BROADENING
            else:
                warnings.warn(f"No broadening data for {diluent} in {molec}. Using default.")
                self.broadening_coefficients[diluent] = self.DEFAULT_BROADENING

    def _download_broadening_file(self, molec, diluent, local_path):
        """Download broadening file from ExoMol if available"""
        url = f"{EXOMOL_URL}{self.molecule}/{molec}/{molec}__{diluent}.broad"
        try:
            with urlopen(url) as response:
                with open(local_path / f"{molec}__{diluent}.broad", 'wb') as f:
                    f.write(response.read())
        except HTTPError as e:
            raise ValueError(f"Broadening file for {diluent} not found at {url}: {e}")

    def get_broadening(self, J, T, Tref=296.0):
        """Compute total broadening for given J and temperature"""
        gamma_total = 0.0
        for diluent, fraction in self.broadening_partners.items():
            coeffs = self.broadening_coefficients.get(diluent, self.DEFAULT_BROADENING)
            gamma = coeffs.get(J, self.DEFAULT_BROADENING)['gamma']
            n = coeffs.get(J, self.DEFAULT_BROADENING)['n']
            gamma_total += fraction * gamma * (Tref / T) ** n
        return gamma_total

    # Existing methods (assumed to be present)
    def _init_file_paths(self, molec, t0):
        pass  # Placeholder for your existing method

    def _load_def_and_pf_files(self, molec):
        pass  # Placeholder for your existing method

    def _load_states_and_transitions(self, molec, skip_optional_data):
        pass  # Placeholder for your existing method

    def QT_interp(self, T):
        """Interpolate partition function at temperature T"""
        return np.interp(T, self.T_gQT, self.gQT)

    def qr_interp(self, T):
        """Partition function ratio Q(T)/Q(Tref)"""
        return self.QT_interp(T) / self.QT_interp(self.Tref)

    def to_partition_function_tabulator(self):
        """Generate PartFuncExoMol object"""
        from radis.levels.partfunc import PartFuncExoMol
        return PartFuncExoMol(self.isotope_fullname, self.T_gQT, self.gQT)
