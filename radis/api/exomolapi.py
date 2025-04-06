"""Molecular database (MDB) class with multi-species broadening support

   * Enhanced MdbExomol handles multiple broadening partners

Key Changes:
- Added broadening_partners parameter to specify multiple species
- Implemented weighted broadening coefficient calculation
- Improved error handling for missing broadening files
- Maintained backward compatibility
"""

import os
import pathlib
import urllib.request
import warnings
import numpy as np
from radis.api.dbmanager import DatabaseManager, get_auto_MEMORY_MAPPING_ENGINE
from radis.db.classes import EXOMOL_MOLECULES, EXOMOL_ONLY_ISOTOPES_NAMES

EXOMOL_URL = "https://www.exomol.com/db/"

# ... [Keep all existing imports and utility functions unchanged] ...

class MdbExomol(DatabaseManager):
    """Molecular database of ExoMol with multi-species broadening support"""
    
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
        broadening_partners=None,  # New: list of broadening partners
    ):
        super().__init__(
            name,
            molecule,
            local_databases,
            engine,
            verbose=verbose,
        )
        assert cache  # Required for caching
        
        # Initialize paths and basic parameters
        self.path = pathlib.Path(path)
        if local_databases is not None:
            self.path = pathlib.Path(local_databases).expanduser() / self.path

        t0 = self.path.parents[0].stem
        molec = t0 + "__" + str(self.path.stem)

        # Store calculation parameters
        self.crit = crit
        self.margin = margin
        self.nurange = [np.min(nurange), np.max(nurange)]
        self.wmin, self.wmax = np.min(nurange), np.max(nurange)
        self.broadf = broadf
        self.broadf_download = broadf_download
        self.engine = engine
        self.broadening_coefficients = {}  # Stores coefficients for all partners

        # New: Handle broadening partners
        self.broadening_partners = self._validate_broadening_partners(
            broadening_partners if broadening_partners is not None else ["air"]
        )

        # Initialize file paths
        self._init_file_paths(molec, t0)
        
        # Load definition and partition function files
        self._load_def_and_pf_files(molec)
        
        # Load states and transition files
        self._load_states_and_transitions(molec, skip_optional_data)

    def _init_file_paths(self, molec, t0):
        """Initialize all required file paths"""
        self.states_file = self.path / pathlib.Path(molec + ".states.bz2")
        self.pf_file = self.path / pathlib.Path(molec + ".pf")
        self.def_file = self.path / pathlib.Path(molec + ".def")
        
        # Broadening files for all supported partners
        self.broad_files = {
            partner: self.path / pathlib.Path(t0 + "__" + partner + ".broad")
            for partner in self.broadening_partners + ["air", "self"]  # Always include these
        }

    def _load_def_and_pf_files(self, molec):
        """Load definition and partition function files"""
        mgr = self.get_datafile_manager()
        
        # Download missing files
        if not self.def_file.exists():
            self.download(molec, extension=[".def"])
        if not self.pf_file.exists():
            self.download(molec, extension=[".pf"])
        
        # Load definition file
        dic_def = read_def(self.def_file)
        self.n_Texp_def = dic_def["n_Texp"] or 0.5  # Default if None
        self.alpha_ref_def = dic_def["alpha_ref"] or 0.07  # Default if None
        self.molmass = dic_def["molmass"]
        
        # Load partition function
        pf = read_pf(self.pf_file)
        self.gQT = pf["QT"].to_numpy()
        self.T_gQT = pf["T"].to_numpy()
        self.Tref = 296.0
        self.QTref = np.array(self.QT_interp(self.Tref))

    def _load_states_and_transitions(self, molec, skip_optional_data):
        """Load states and transition data"""
        mgr = self.get_datafile_manager()
        
        # Download and cache states file
        if not self.states_file.exists() and not mgr.cache_file(self.states_file).exists():
            self.download(molec, extension=[".states.bz2"])
        
        # Load states data
        states = self._load_cached_file(
            self.states_file,
            lambda: read_states(
                self.states_file,
                read_def(self.def_file),
                engine="vaex" if self.engine == "vaex" else "csv",
                skip_optional_data=skip_optional_data,
            )
        )

        # Handle transition files
        dic_def = read_def(self.def_file)
        if dic_def["numinf"] is None:
            self.trans_file = [self.path / pathlib.Path(molec + ".trans.bz2")]
            self.num_tag = [None]
        else:
            self._init_multipart_transition_files(molec, dic_def)

        # Download broadening files if needed
        self._download_broadening_files(molec)

        # Load and cache transition files
        for trans_file, num_tag in zip(self.trans_file, self.num_tag):
            self._load_and_cache_transitions(trans_file, num_tag, states, dic_def, skip_optional_data)

    def _load_and_cache_transitions(self, trans_file, num_tag, states, dic_def, skip_optional_data):
        """Load and cache transition data"""
        mgr = self.get_datafile_manager()
        
        if not trans_file.exists():
            self.download(molec, extension=[".trans.bz2"], numtag=num_tag)
        
        trans = self._load_cached_file(
            trans_file,
            lambda: read_trans(trans_file, engine="vaex" if self.engine == "vaex" else "csv"),
            postprocess=lambda df: pickup_gE(states, df, dic_def, skip_optional_data, self.engine)
        )
        
        # Calculate line strengths
        from radis.lbl.base import linestrength_from_Einstein
        trans["Sij0"] = linestrength_from_Einstein(
            A=trans["A"],
            gu=trans["gup"],
            El=trans["elower"],
            Ia=1,
            nu=trans["nu_lines"],
            Q=self.QTref,
            T=self.Tref,
        )

    def set_broadening_coef(
        self,
        df,
        alpha_ref_def=None,
        n_Texp_def=None,
        output=None,
        add_columns=True,
        species=None,
        partner_mole_fractions=None,
    ):
        """Calculate effective broadening from multiple partners
        
        Parameters
        ----------
        species : str or list
            Broadening partner(s) to use
        partner_mole_fractions : dict
            Dictionary of {partner: mole_fraction}
        """
        if species is None:
            species = self.broadening_partners
            
        if isinstance(species, str):
            species = [species]
            
        # Validate partners
        species = self._validate_broadening_partners(species)
        
        # Calculate mole fractions
        if partner_mole_fractions is None:
            partner_mole_fractions = {s: 1/len(species) for s in species}
        elif abs(sum(partner_mole_fractions.values()) - 1) > 1e-6:
            raise ValueError("Mole fractions must sum to 1")
        
        # Calculate effective broadening parameters
        effective_alpha_ref = 0
        effective_n_Texp = 0
        
        for partner in species:
            alpha_ref, n_Texp = self._get_broadening_coefficients(df, partner)
            effective_alpha_ref += partner_mole_fractions[partner] * alpha_ref
            effective_n_Texp += partner_mole_fractions[partner] * n_Texp
        
        # Store results
        self.alpha_ref = effective_alpha_ref
        self.n_Texp = effective_n_Texp
        
        if add_columns:
            self._add_broadening_columns(df, species, partner_mole_fractions)

    def _get_broadening_coefficients(self, df, partner):
        """Get broadening coefficients for a specific partner"""
        if partner not in self.broadening_coefficients:
            self._fetch_broadening_coefficients(df, partner)
        return self.broadening_coefficients[partner]

    def _fetch_broadening_coefficients(self, df, partner):
        """Fetch and calculate broadening coefficients for a partner"""
        file = self.broad_files.get(partner)
        
        # Use defaults for air/self if file not found
        if not file or not os.path.exists(file):
            if partner in ["air", "self"]:
                self.broadening_coefficients[partner] = (
                    np.full(len(df), self.alpha_ref_def),
                    np.full(len(df), self.n_Texp_def)
                )
                return
            raise ValueError(f"Broadening file for {partner} not found")
        
        bdat = read_broad(file, output=self.engine)
        codelv = check_code_level(bdat, output=self.engine)
        
        if codelv == "a0":
            alpha_ref, n_Texp = make_j2b(
                bdat,
                alpha_ref_default=self.alpha_ref_def,
                n_Texp_default=self.n_Texp_def,
                jlower_max=df["jlower"].max(),
                output=self.engine,
            )
        elif codelv == "a1":
            j2alpha_ref, j2n_Texp = make_j2b(
                bdat,
                alpha_ref_default=self.alpha_ref_def,
                n_Texp_default=self.n_Texp_def,
                jlower_max=df["jlower"].max(),
                output=self.engine,
            )
            alpha_ref, n_Texp = make_jj2b(
                bdat,
                j2alpha_ref_def=j2alpha_ref,
                j2n_Texp_def=j2n_Texp,
                jupper_max=df["jupper"].max(),
                output=self.engine,
            )
        else:
            alpha_ref = np.full(len(df), self.alpha_ref_def)
            n_Texp = np.full(len(df), self.n_Texp_def)
            
        self.broadening_coefficients[partner] = (alpha_ref, n_Texp)

    def _add_broadening_columns(self, df, species, mole_fractions):
        """Add broadening columns to DataFrame based on partners"""
        if len(species) == 1 and species[0] == "air":
            # Backward compatibility
            self.add_column(df, "airbrd", self.alpha_ref)
            self.add_column(df, "Tdpair", self.n_Texp)
        else:
            # New format for multiple partners
            self.add_column(df, "effective_brd", self.alpha_ref)
            self.add_column(df, "effective_Tdp", self.n_Texp)
            # Store partner info as metadata
            df.attrs["broadening_partners"] = species
            df.attrs["broadening_mole_fractions"] = mole_fractions

    def _validate_broadening_partners(self, partners):
        """Check which broadening partners are available"""
        available = []
        for partner in partners:
            if partner in ["air", "self"]:
                available.append(partner)
            elif (self.path / pathlib.Path(self.isotope_fullname + "__" + partner + ".broad")).exists():
                available.append(partner)
            elif self.verbose:
                print(f"Warning: Broadening coefficients for {partner} not found")
        
        if not available:
            raise ValueError("No valid broadening partners found")
        return available

    def _download_broadening_files(self, molec):
        """Download required broadening files"""
        if not self.broadf or not self.broadf_download:
            return
            
        # Check if we need to download any files
        need_download = any(
            not f.exists() and p not in ["air", "self"]
            for p, f in self.broad_files.items()
        )
        
        if need_download:
            self.download(molec, extension=[".broad"])

    # ... [Keep all other existing methods unchanged] ...

# ... [Keep all remaining utility functions unchanged] ...
