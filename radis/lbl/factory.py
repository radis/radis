# -*- coding: utf-8 -*-
"""Contains the :py:class:`~radis.lbl.factory.SpectrumFactory` class, which is
the core of the RADIS Line-by-Line module.
[Existing docstring remains unchanged]
"""

from typing import Union, Dict
from warnings import warn

import astropy.units as u
import numpy as np
from numpy import arange, exp, expm1
from scipy.optimize import OptimizeResult

from radis import version
from radis.db import MOLECULES_LIST_EQUILIBRIUM, MOLECULES_LIST_NONEQUILIBRIUM
from radis.db.classes import get_molecule, get_molecule_identifier, to_conventional_name

try:
    from .bands import BandFactory
    from .base import get_wavenumber_range
except ImportError:
    from radis.lbl.bands import BandFactory
    from radis.lbl.base import get_wavenumber_range

from radis import config
from radis.db.classes import is_atom, is_neutral
from radis.misc.basics import flatten, is_float, is_range, list_if_float, round_off
from radis.misc.utils import Default
from radis.phys.constants import k_b
from radis.phys.convert import conv2
from radis.phys.units import convert_universal
from radis.phys.units_astropy import convert_and_strip_units
from radis.spectrum.equations import calc_radiance
from radis.spectrum.spectrum import Spectrum

class SpectrumFactory(BandFactory):
    """A class to put together all functions related to loading CDSD / HITRAN
    databases, calculating the broadenings, and summing over all the lines.
    [Existing docstring remains unchanged]
    """

    __slots__ = BandFactory.__slots__

    def __init__(
        self,
        wmin=None,
        wmax=None,
        wunit=Default("cm-1"),
        wavenum_min=None,
        wavenum_max=None,
        wavelength_min=None,
        wavelength_max=None,
        Tref=296,
        pressure=1.01325,
        mole_fraction=1,
        path_length=1,
        wstep=0.01,
        isotope="all",
        medium="air",
        truncation=Default(50),
        neighbour_lines=0,
        pseudo_continuum_threshold=0,
        self_absorption=True,
        chunksize=None,
        optimization="simple",
        folding_thresh=1e-6,
        zero_padding=-1,
        broadening_method="voigt",
        cutoff=0,
        parsum_mode="full summation",
        verbose=True,
        warnings=True,
        species=None,
        save_memory=False,
        export_populations=None,
        export_lines=False,
        gpu_backend=None,
        diluent=Default(None),
        broadening_partners=None,  # NEW: Add broadening_partners parameter
        lbfunc=None,
        potential_lowering=None,
        pfsource="default",
        **kwargs,
    ):
        # Initialize BandFactory object
        super(SpectrumFactory, self).__init__()

        # Check inputs (deal with deprecated format)
        if medium not in ["air", "vacuum"]:
            raise ValueError("Wavelength must be one of: 'air', 'vacuum'")
        kwargs0 = kwargs
        if "molecule" in kwargs:
            if species is not None and species != kwargs["molecule"]:
                raise Exception(
                    "Both `molecule` and `species` arguments have been given and aren't equal, but `molecule` is just an alias for `species`."
                )
            species = molecule = kwargs["molecule"]
            kwargs0.pop("molecule")

        # [Existing input checks remain unchanged]

        # Handle broadening_partners and diluent compatibility
        if broadening_partners is not None:
            if not isinstance(broadening_partners, dict):
                raise ValueError("broadening_partners must be a dictionary of {diluent: mole_fraction}")
            if diluent is not Default(None):
                warn("Both `diluent` and `broadening_partners` specified. Using `broadening_partners`.")
            diluent = broadening_partners  # Use broadening_partners as the diluent spec
        elif isinstance(diluent, Default):
            # Default diluent logic remains if broadening_partners not provided
            diluent = None  # Will be set later based on molecule type

        # Calculate waveranges
        wavenum_min, wavenum_max, input_wunit = get_wavenumber_range(
            wmin, wmax, wunit, wavenum_min, wavenum_max, wavelength_min, wavelength_max, medium, return_input_wunit=True
        )
        self.input_wunit = input_wunit
        self._wstep = wstep

        # Set default variables from config
        self._sparse_ldm = config["SPARSE_WAVERANGE"]
        self.params["sparse_ldm"] = config["SPARSE_WAVERANGE"]

        # Init variables
        if molecule is not None:
            species = molecule = to_conventional_name(molecule)
            if isinstance(molecule, int):
                species = molecule = get_molecule(molecule)
            if not is_atom(molecule):
                if molecule not in MOLECULES_LIST_EQUILIBRIUM + MOLECULES_LIST_NONEQUILIBRIUM:
                    raise ValueError(f"Unsupported molecule: {molecule}")
                self.input.isatom = False
                self.input.isneutral = None
                if diluent is None:  # Only set default if not overridden by broadening_partners
                    diluent = "air"
            else:
                self.input.isatom = True
                self.input.isneutral = is_neutral(molecule)
                if diluent is None:  # Only set default if not overridden by broadening_partners
                    diluent = "H"

        # Store isotope identifier
        if not isinstance(isotope, str):
            isotope = ",".join([str(k) for k in list_if_float(isotope)])

        # Check molecule in diluent
        if isinstance(diluent, str) and diluent == molecule:
            raise KeyError(f"{molecule} cannot be both molecule and diluent.")
        elif isinstance(diluent, dict) and molecule in diluent:
            raise KeyError(f"{molecule} cannot be both molecule and a diluent in broadening_partners.")

        # Initialize input conditions
        self.input.wavenum_min = wavenum_min
        self.input.wavenum_max = wavenum_max
        self.input.Tref = convert_and_strip_units(Tref, u.K)
        self.input.pressure = convert_and_strip_units(pressure, u.bar)
        self.input.mole_fraction = mole_fraction
        self.dataframe_engine = config["DATAFRAME_ENGINE"] if config["DATAFRAME_ENGINE"] in ["pandas", "vaex"] else "pandas"
        self.input.path_length = convert_and_strip_units(path_length, u.cm)
        self.input.species = species
        self.input.state = "X"
        self.input.isotope = isotope
        self.input.self_absorption = self_absorption
        self.input.potential_lowering = potential_lowering
        self.input.pfsource = pfsource

        # Initialize computation variables
        self.params.wstep = wstep
        self.params.pseudo_continuum_threshold = pseudo_continuum_threshold
        self.params.diluent = diluent  # Now can be a dict from broadening_partners or str
        self.params.broadening_partners = broadening_partners  # Store explicitly for reference
        self._diluent = None  # Computed later in _generate_diluent_molefraction

        # [Remaining initialization logic unchanged]
        if cutoff is None:
            cutoff = 0
        self.params.cutoff = cutoff
        self.params.parsum_mode = parsum_mode
        self.verbose = verbose
        self.params.truncation = self.truncation = truncation.value if isinstance(truncation, Default) else truncation
        self.params.neighbour_lines = neighbour_lines
        self.params.wavenum_min_calc = wavenum_min - neighbour_lines
        self.params.wavenum_max_calc = wavenum_max + neighbour_lines
        self._neighbour_lines = neighbour_lines
        self.params.broadening_method = broadening_method
        self.params.optimization = optimization
        self.params.folding_thresh = folding_thresh
        self.misc.zero_padding = zero_padding
        self.params.lbfunc = lbfunc
        self.misc.chunksize = chunksize
        self.save_memory = save_memory
        self.autoupdatedatabase = False
        self.autoretrievedatabase = False
        self.SpecDatabase = None
        self.database = None
        if isinstance(warnings, dict):
            self.warnings.update(warnings)
        elif warnings in [True, "warn", "warning"]:
            self.warnings["default"] = "warn"
        elif warnings == "error":
            self.warnings["default"] = "error"
        elif warnings in [False, "ignore"]:
            self.warnings["default"] = "ignore"
        else:
            raise ValueError(f"Unexpected value for warnings: {warnings}")
        self.misc.warning_linestrength_cutoff = 1e-2
        self.misc.warning_broadening_threshold = 1e-2

    def eq_spectrum(
        self,
        Tgas,
        mole_fraction=None,
        path_length=None,
        diluent=None,
        broadening_partners=None,  # NEW: Add broadening_partners as optional override
        pressure=None,
        name=None,
    ) -> Spectrum:
        """Generate a spectrum at equilibrium."""
        # [Existing preprocessing unchanged until diluent handling]
        Tgas = convert_and_strip_units(Tgas, u.K)
        path_length = convert_and_strip_units(path_length, u.cm)
        pressure = convert_and_strip_units(pressure, u.bar)
        if path_length is not None:
            self.input.path_length = path_length
        if mole_fraction is not None:
            self.input.mole_fraction = mole_fraction
        if pressure is not None:
            self.input.pressure = pressure
        if not is_float(Tgas):
            raise ValueError(f"Tgas should be float or Astropy unit. Got {Tgas}")
        self.input.Tgas = Tgas
        pressure = self.input.pressure
        mole_fraction = self.input.mole_fraction
        path_length = self.input.path_length
        verbose = self.verbose
        self._reset_profiler(verbose)
        self._check_inputs(mole_fraction, max(flatten(Tgas)))
        if self.autoretrievedatabase:
            s = self._retrieve_from_database()
            if s is not None:
                return s

        self.profiler.start("spectrum_calculation", 1)
        self.profiler.start("spectrum_calc_before_obj", 2)
        self._check_line_databank()
        self._reinitialize()

        self.calc_linestrength_eq(Tgas)
        self._cutoff_linestrength()

        # Handle broadening_partners override
        if broadening_partners is not None:
            if not isinstance(broadening_partners, dict):
                raise ValueError("broadening_partners must be a dictionary of {diluent: mole_fraction}")
            diluent = broadening_partners
        self._generate_diluent_molefraction(mole_fraction, diluent)

        self._calc_broadening_HWHM()
        self.calc_lineshift()
        self._generate_wavenumber_arrays()
        I_continuum = self.calculate_pseudo_continuum()
        wavenumber, abscoeff_v = self._calc_broadening()
        abscoeff_v = self._add_pseudo_continuum(abscoeff_v, I_continuum)

        # [Remaining calculation and export logic unchanged]
        density = mole_fraction * ((pressure * 1e5) / (k_b * Tgas)) * 1e-6
        abscoeff = abscoeff_v * density
        absorbance = abscoeff * path_length
        emissivity_noslit = -expm1(-absorbance)
        transmittance_noslit = 1 - emissivity_noslit
        radiance_noslit = calc_radiance(wavenumber, emissivity_noslit, Tgas, unit=self.units["radiance_noslit"])
        assert self.units["abscoeff"] == "cm-1"
        # [Export logic remains unchanged]
        conditions = self.get_conditions(add_config=True)
        conditions.update({
            "calculation_time": self.profiler.final[list(self.profiler.final)[-1]]["spectrum_calc_before_obj"],
            "lines_calculated": self._Nlines_calculated,
            "lines_cutoff": self._Nlines_cutoff,
            "lines_in_continuum": self._Nlines_in_continuum,
            "thermal_equilibrium": True,
            "diluents": self._diluent,
            "radis_version": version,
            "spectral_points": int((self.params.wavenum_max_calc - self.params.wavenum_min_calc) / self.params.wstep),
            "profiler": dict(self.profiler.final),
        })
        if self.params.optimization is not None:
            conditions.update({"NwL": self.NwL, "NwG": self.NwG})
        del self.profiler.final[list(self.profiler.final)[-1]]["spectrum_calc_before_obj"]
        populations = None
        lines = self.get_lines()
        quantities = {
            "wavenumber": wavenumber,
            "abscoeff": abscoeff,
            "absorbance": absorbance,
            "emissivity_noslit": emissivity_noslit,
            "transmittance_noslit": transmittance_noslit,
            "radiance_noslit": radiance_noslit,
        }
        if I_continuum is not None and self._export_continuum:
            quantities.update({"abscoeff_continuum": I_continuum * density})
        conditions["default_output_unit"] = self.input_wunit
        s = Spectrum(
            quantities=quantities,
            conditions=conditions,
            populations=populations,
            lines=lines,
            units=self.units,
            cond_units=self.cond_units,
            check_wavespace=False,
            name=name,
            references=dict(self.reftracker),
        )
        if self.autoupdatedatabase:
            self.SpecDatabase.add(s, if_exists_then="increment")
        self.profiler.stop("generate_spectrum_obj", "Generated Spectrum object")
        self.profiler.stop("spectrum_calculation", "Spectrum calculated")
        return s

    def _generate_diluent_molefraction(self, mole_fraction, diluent):
        """Generate diluent mole fractions based on input diluent or broadening_partners."""
        from radis.misc.warning import MoleFractionError

        molecule = self.input.species
        if diluent is None:
            # Use factory default if not overridden
            if isinstance(self.params.diluent, str):
                diluents = {self.params.diluent: 1 - mole_fraction}
            elif isinstance(self.params.diluent, dict):
                diluents = self.params.diluent.copy()
            else:
                diluents = {"air": 1 - mole_fraction}  # Fallback default
        else:
            if isinstance(diluent, str):
                diluents = {diluent: 1 - mole_fraction}
            elif isinstance(diluent, dict):
                diluents = diluent.copy()
            else:
                raise ValueError("diluent must be a string or dictionary")

        # Check for molecule in diluents
        if molecule in diluents:
            raise KeyError(f"{molecule} cannot be both molecule and a diluent.")

        # Validate total mole fraction
        total_mole_fraction = mole_fraction + sum(diluents.values())
        if not np.isclose(total_mole_fraction, 1, rtol=1e-5):
            if total_mole_fraction > 1:
                message = f"Total mole fraction = {total_mole_fraction} exceeds 1. Adjust mole_fraction or diluents."
            else:
                message = f"Total mole fraction = {total_mole_fraction} less than 1. Adjust mole_fraction or diluents."
            raise MoleFractionError(message)
        self._diluent = diluents

    # [Remaining methods like non_eq_spectrum, optically_thin_power, etc., would need similar updates
    #  to accept broadening_partners and pass it to _generate_diluent_molefraction, but are omitted
    #  for brevity. Only key changes are shown above.]

# [Rest of the file (extra functions, test section) remains unchanged]
