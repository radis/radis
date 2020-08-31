# -*- coding: utf-8 -*-
"""
Summary
-------

Public (front-end) functions to calculate Spectrum with HITRAN / CDSD databanks.
Uses the SpectrumFactory classe from `factory.py`, Spectrum from `spectrum.py`
and line survey from `line_survey.py` 

Routine Listing
---------------

:func:`~radis.lbl.calc.calc_spectrum`

-------------------------------------------------------------------------------

"""

from __future__ import print_function, absolute_import, division, unicode_literals

from copy import deepcopy
from radis.lbl.factory import SpectrumFactory
from radis.phys.convert import nm2cm
from radis.misc.basics import all_in
from os.path import exists

# %%
def calc_spectrum(
    wavenum_min=None,
    wavenum_max=None,
    wavelength_min=None,
    wavelength_max=None,
    Tgas=None,
    Tvib=None,
    Trot=None,
    pressure=1.01325,
    molecule=None,
    isotope="all",
    mole_fraction=1,
    path_length=1,
    medium="air",
    databank="fetch",
    wstep=0.01,
    broadening_max_width=10,
    optimization="min-RMS",
    overpopulation=None,
    name=None,
    use_cached=True,
    verbose=True,
    mode="cpu",
    **kwargs
):
    """Multipurpose function to calculate :class:`~radis.spectrum.spectrum.Spectrum`
        under equilibrium (using either CPU or GPU), or non-equilibrium, with or without overpopulation.
        It's a wrapper to :class:`~radis.lbl.factory.SpectrumFactory` class.
        For advanced used, please refer to the aforementionned class.
    ​
        Parameters
        ----------
        wavenum_min: float [cm-1]
            minimum wavenumber to be processed in cm^-1
        wavenum_max: float [cm-1]
            maximum wavenumber to be processed in cm^-1
    ​
        wavelength_min: float [nm]
            minimum wavelength to be processed in nm. Wavelength in ``'air'`` or
            ``'vacuum'`` depending of the value of the parameter ``'medium='``
    ​
        wavelength_max: float [nm]
            maximum wavelength to be processed in nm. Wavelength in ``'air'`` or
            ``'vacuum'`` depending of the value of the parameter ``'medium='``
    ​
        Tgas: float [K]
            Gas temperature. If non equilibrium, is used for Ttranslational.
            Default ``300`` K
    ​
        Tvib: float [K]
            Vibrational temperature. If ``None``, equilibrium calculation is run with Tgas
    ​
        Trot: float [K]
            Rotational temperature. If ``None``, equilibrium calculation is run with Tgas
    ​
        pressure: float [bar]
            partial pressure of gas in bar. Default ``1.01325`` (1 atm)
    ​
        molecule: int, str, list or ``None``
            molecule id (HITRAN format) or name. For multiple molecules, use a list.
            The `isotope`, `mole_fraction`, `databank` and `overpopulation` parameters must then
            be dictionaries.
            If ``None``, the molecule can be infered
            from the database files being loaded. See the list of supported molecules
            in :py:data:`~radis.io.MOLECULES_LIST_EQUILIBRIUM`
            and :py:data:`~radis.io.MOLECULES_LIST_NONEQUILIBRIUM`.
            Default ``None``.
    ​
        isotope: int, list, str of the form ``'1,2'``, or ``'all'``, or dict
            isotope id (sorted by relative density: (eg: 1: CO2-626, 2: CO2-636 for CO2).
            See [HITRAN-2016]_ documentation for isotope list for all species. If ``'all'``,
            all isotopes in database are used (this may result in larger computation
            times!). Default ``'all'``.

            For multiple molecules, use a dictionary with molecule names as keys.

            Example::

                mole_fraction={'CO2':0.8 ,  'CO':0.2 }
    ​
        mole_fraction: float or dict
            database species mole fraction. Default ``1``.

            For multiple molecules, use a dictionary with molecule names as keys.

            Example::

                mole_fraction={'CO2': 0.8, 'CO':0.2}
    ​
        path_length: float [cm]
            slab size. Default ``1``.
    ​
        databank: str or dict
            can be either:
    ​
            - ``'fetch'``, to fetch automatically from [HITRAN-2016]_ through astroquery.

            .. warning::

                [HITRAN-2016]_ is valid for low temperatures (typically < 700 K). For higher
                temperatures you may need [HITEMP-2010]_

            - the name of a a valid database file, in which case the format is inferred.
              For instance, ``'.par'`` is recognized as ``hitran/hitemp`` format.
              Accepts wildcards ``'*'`` to select multiple files.

            - the name of a spectral database registered in your ``~/.radis``
              configuration file. This allows to use multiple database files.
              See :ref:`Configuration file <label_lbl_config_file>`.
    ​
            Default ``'fetch'``. See :class:`~radis.lbl.loader.DatabankLoader` for more
            information on line databases, and :data:`~radis.misc.config.DBFORMAT` for
            your ``~/.radis`` file format


            For multiple molecules, use a dictionary with molecule names as keys::
    ​
                databank='fetch'     # automatic download
                databank='PATH/TO/05_HITEMP2019.par'    # path to a file
                databank='*CO2*.par' #to get all the files that have CO2 in their names (case insensitive)
                databank='HITEMP-2019-CO'   # user-defined database in Configuration file
                databank = {'CO2' : 'PATH/TO/05_HITEMP2019.par', 'CO' : 'fetch'}  # for multiple molecules
    ​
        medium: ``'air'``, ``'vacuum'``
            propagating medium when giving inputs with ``'wavenum_min'``, ``'wavenum_max'``.
            Does not change anything when giving inputs in wavenumber. Default ``'air'``
    ​
        wstep: float (cm-1)
            Spacing of calculated spectrum. Default ``0.01 cm-1``
    ​
        broadening_max_width: float (cm-1)
            Full width over which to compute the broadening. Large values will create
            a huge performance drop (scales as ~broadening_width^2 without DLM)
            The calculated spectral range is increased (by broadening_max_width/2
            on each side) to take into account overlaps from out-of-range lines.
            Default ``10`` cm-1.
    ​
        Other Parameters
        ----------------
    ​
        optimization : ``"simple"``, ``"min-RMS"``, ``None``
            If either ``"simple"`` or ``"min-RMS"`` DLM optimization for lineshape calculation is used:
            - ``"min-RMS"`` : weights optimized by analytical minimization of the RMS-error (See: [DLM_article]_)
            - ``"simple"`` : weights equal to their relative position in the grid

            If using the DLM optimization, broadening method is automatically set to ``'fft'``.
            If ``None``, no lineshape interpolation is performed and the lineshape of all lines is calculated.

            Refer to [DLM_article]_ for more explanation on the DLM method for lineshape interpolation.

            Default ``"min-RMS"``
    ​
        overpopulation: dict
            dictionary of overpopulation compared to the given vibrational temperature.
            Default ``None``.

            Example::

                overpopulation = {'CO2' : {'(00`0`0)->(00`0`1)': 2.5,
                                           '(00`0`1)->(00`0`2)': 1,
                                           '(01`1`0)->(01`1`1)': 1,
                                           '(01`1`1)->(01`1`2)': 1

                                            }
                                 }
    ​
        slit: float, str, or ``None``
            if float, FWHM of a triangular slit function. If str, path to an
            experimental slit function. If None, no slit is applied. Default ``None``.
    ​
        plot: str
            any parameter such as 'radiance' (if slit is given), 'radiance_noslit',
            'absorbance', etc...   Default ``None``
    ​
        name: str
            name of the case. If None, a unique ID is generated. Default ``None``
    ​
        use_cached: boolean
            use cached files for line database and energy database. Default ``True``
    ​
        verbose: boolean, or int
            If ``False``, stays quiet. If ``True``, tells what is going on.
            If ``>=2``, gives more detailed messages (for instance, details of
            calculation times). Default ``True``.
    ​
        **kwargs: other inputs forwarded to SpectrumFactory
            For instance: ``warnings``.
            See :class:`~radis.lbl.factory.SpectrumFactory` documentation for more
            details on input.
            For instance:
    ​
        pseudo_continuum_threshold: float
            if not 0, first calculate a rough approximation of the spectrum, then
            moves all lines whose linestrength intensity is less than this threshold
            of the maximum in a semi-continuum. Values above 0.01 can yield significant
            errors, mostly in highly populated areas. 80% of the lines can typically
            be moved in a continuum, resulting in 5 times faster spectra. If 0,
            no semi-continuum is used. Default 0.

        mode: ``'cpu'``, ``'gpu'``
            if set to 'cpu', computes the spectra purely on the CPU. if set to 'gpu',
            offloads the calculations of lineshape and broadening steps to the GPU
            making use of parallel computations to speed up the process. Default 'cpu'.
            Note that mode='gpu' requires CUDA compatible hardware to execute. For more information on how to setup your system to run GPU-accelerated methods using CUDA and Cython, check `GPU Spectrum Calculation on RADIS <https://radis.readthedocs.io/en/latest/lbl/gpu.html>`
    ​
        Returns
        -------

        s: :class:`~radis.spectrum.spectrum.Spectrum`
            Output spectrum.
    ​
            Use the :py:meth:`~radis.spectrum.spectrum.Spectrum.get` method to retrieve a
            spectral quantity (``'radiance'``, ``'radiance_noslit'``, ``'absorbance'``, etc...)
    ​
            Or the :py:meth:`~radis.spectrum.spectrum.Spectrum.plot` method to plot it
            directly.
    ​
            See [1]_ to get an overview of all Spectrum methods
    ​
        References
        ----------
    ​
        .. [1] RADIS doc: `Spectrum how to? <https://radis.readthedocs.io/en/latest/spectrum/spectrum.html#label-spectrum>`__
    ​    .. [2] RADIS GPU support: 'GPU Calculations on RADIS <https://radis.readthedocs.io/en/latest/lbl/gpu.html>'
    ​
        Examples
        --------
    ​
        Calculate a CO spectrum from the HITRAN database::
    ​
            s = calc_spectrum(1900, 2300,         # cm-1
                              molecule='CO',
                              isotope='1,2,3',
                              pressure=1.01325,   # bar
                              Tgas=1000,
                              mole_fraction=0.1,
                              )
            s.apply_slit(0.5, 'nm')
            s.plot('radiance')

        This example uses the :py:meth:`~radis.spectrum.spectrum.Spectrum.apply_slit`
        and :py:meth:`~radis.spectrum.spectrum.Spectrum.plot` methods. See also
        :py:meth:`~radis.spectrum.spectrum.Spectrum.line_survey`::

            s.line_survey(overlay='radiance')

        Calculate a CO2 spectrum from the CDSD-4000 database:

            s = calc_spectrum(2200, 2400,   # cm-1
                              molecule='CO2',
                              isotope='1',
                              databank='/path/to/cdsd/databank/in/npy/format/',
                              pressure=0.1,  # bar
                              Tgas=1000,
                              mole_fraction=0.1,
                              mode='gpu'
                              )

            s.plot('absorbance')

        This example uses the :py:meth:`~radis.lbl.factor.eq_spectrum_gpu` method to calculate
        the spectrum on the GPU. The databank points to the CDSD-4000 databank that has been
        pre-processed and stored in `numpy.npy` format.
    ​
        Refer to the online :ref:`Examples <label_examples>` for more cases, and to
        the :ref:`Spectrum page <label_spectrum>` for details on post-processing methods.

        For more details on how to use the GPU method and process the database, refer to the examples
        linked above and the documentation on :ref:`GPU support for RADIS <label_gpu>`.
    ​
        See Also
        --------

        :class:`~radis.lbl.factory.SpectrumFactory`,
        the :ref:`Spectrum page <label_spectrum>`
    """

    from radis.los.slabs import MergeSlabs

    if molecule is not None and type(molecule) != list:
        molecule = [molecule]  # fall back to the other case: multiple molecules

    # Stage 1. Find all molecules, whatever the user input configuration

    # ... Input arguments that CAN be dictionaries of molecules.
    DICT_INPUT_ARGUMENTS = {
        "isotope": isotope,
        "mole_fraction": mole_fraction,
        "databank": databank,
    }
    # Same, but when the values of the arguments themselves are already a dict.
    # (dealt with separately because we cannot use them to guess what are the input molecules)
    DICT_INPUT_DICT_ARGUMENTS = {"overpopulation": overpopulation}

    def _check_molecules_are_consistent(
        molecule_reference_set, reference_name, new_argument, new_argument_name
    ):
        """Will test that molecules set are the same in molecule_reference_set and new_argument, if new_argument is a dict.
        molecule_reference_set is a set of molecules (yeah!).
        reference_name is the name of the argument from which we guessed the list of molecules (used to have a clear error message).
        new_argument is the new argument to check
        new_argument_name is its name

        Returns the set of molecules as found in new_argument, if applicable, else the molecule_reference_set (this allows us to parse all arguments too)

        Note that names are just here to provide clear error messages to the user if there is a contradiction.
        """

        if isinstance(new_argument, dict):
            if molecule_reference_set is None:  # input molecules are still unknown
                return set(new_argument.keys()), new_argument_name
            elif set(new_argument.keys()) != set(molecule_reference_set):
                raise ValueError(
                    "Keys of molecules in the {0} dictionary must be the same as given in `{1}=`, i.e: {2}. Instead, we got {3}".format(
                        new_argument_name,
                        reference_name,
                        molecule_reference_set,
                        set(new_argument.keys()),
                    )
                )
            else:
                return (
                    set(new_argument.keys()),
                    new_argument_name,
                )  # so now we changed the reference
        else:
            return molecule_reference_set, reference_name

    # Parse all inputs:
    molecule_reference_set = molecule
    reference_name = "molecule"
    for argument_name, argument in DICT_INPUT_ARGUMENTS.items():
        molecule_reference_set, reference_name = _check_molecules_are_consistent(
            molecule_reference_set, reference_name, argument, argument_name
        )

    # ... Now we are sure there are no contradctions. Just ensure we have molecules:
    if molecule_reference_set is None:
        raise ValueError(
            "Please enter the molecule(s) to calculate in the `molecule=` argument or as a dictionary in the following: {0}".format(
                list(DICT_INPUT_ARGUMENTS.keys())
            )
        )

    # Stage 2. Now we have the list of molecules. Let's get the input arguments for each of them.

    # ... Initialize and fill the master-dictionary
    molecule_dict = {}

    for molecule in molecule_reference_set:
        molecule_dict[molecule] = {}

        for argument_name, argument_dict in DICT_INPUT_ARGUMENTS.items():
            if isinstance(argument_dict, dict):
                # Choose the correspond value
                molecule_dict[molecule][argument_name] = argument_dict[molecule]
                # Will raise a KeyError if not defined. That's fine!
                # TODO: maybe need to catch KeyError and raise a better error message?
            else:  # argument_name is not a dictionary.
                # Let's distribute the same value to every molecule:
                molecule_dict[molecule][argument_name] = argument_dict
                # If wrong argument, it will be caught in _calc_spectrum() later.

    # ... Special case of dictionary arguments. Find out if they were given as default, or per dictionary of molecules:
    is_same_for_all_molecules = dict.fromkeys(DICT_INPUT_DICT_ARGUMENTS)
    for argument_name, argument_dict in DICT_INPUT_DICT_ARGUMENTS.items():
        if not isinstance(argument_dict, dict):
            is_same_for_all_molecules[argument_name] = True
        else:
            # Argument is a dictionary. Guess if keys are molecules, or levels.
            # Ex: overpopulation dict could be {'CO2':{'(0,0,0,1)':10}} or directly {{'(0,0,0,1)':10}}
            argument_keys = set(argument_dict.keys())
            if all_in(argument_keys, molecule_reference_set):
                is_same_for_all_molecules[argument_name] = False
            else:
                is_same_for_all_molecules[argument_name] = True
    # ... now fill them in:
    for argument_name, argument_dict in DICT_INPUT_DICT_ARGUMENTS.items():
        if is_same_for_all_molecules[argument_name]:
            for mol in molecule_reference_set:
                # copy the value for everyone
                molecule_dict[mol][argument_name] = deepcopy(
                    argument_dict
                )  # in case it gets edited.
        else:  # argument_dict keys are the molecules:
            for mol in molecule_reference_set:
                molecule_dict[mol][argument_name] = argument_dict[mol]

    # Stage 3: Now let's calculate all the spectra
    s_list = []
    for molecule, dict_arguments in molecule_dict.items():
        kwargs_molecule = deepcopy(
            kwargs
        )  # these are the default supplementary arguments. Deepcopy ensures that they remain the same for all molecules, even if modified in _calc_spectrum

        # We add all of the DICT_INPUT_ARGUMENTS values:
        kwargs_molecule.update(**dict_arguments)

        s_list.append(
            _calc_spectrum(
                wavenum_min=wavenum_min,
                wavenum_max=wavenum_max,
                wavelength_min=wavelength_min,
                wavelength_max=wavelength_max,
                Tgas=Tgas,
                Tvib=Tvib,
                Trot=Trot,
                pressure=pressure,
                # overpopulation=overpopulation,  # now in dict_arguments
                molecule=molecule,
                # isotope=isotope,                # now in dict_arguments
                # mole_fraction=mole_fraction,    # now in dict_arguments
                path_length=path_length,
                # databank=databank,              # now in dict_arguments
                medium=medium,
                wstep=wstep,
                broadening_max_width=broadening_max_width,
                optimization=optimization,
                name=name,
                use_cached=use_cached,
                verbose=verbose,
                mode=mode,
                **kwargs_molecule
            )
        )

    # Stage 4: merge all molecules and return
    return MergeSlabs(*s_list)


def _calc_spectrum(
    wavenum_min,
    wavenum_max,
    wavelength_min,
    wavelength_max,
    Tgas,
    Tvib,
    Trot,
    pressure,
    overpopulation,
    molecule,
    isotope,
    mole_fraction,
    path_length,
    databank,
    medium,
    wstep,
    broadening_max_width,
    optimization,
    name,
    use_cached,
    verbose,
    mode,
    **kwargs
):
    """See :py:func:`~radis.lbl.calc.calc_spectrum`"""

    # Check inputs

    # ... wavelengths / wavenumbers
    if (wavelength_min is not None or wavelength_max is not None) and (
        wavenum_min is not None or wavenum_max is not None
    ):
        raise ValueError("Wavenumber and Wavelength both given... it's time to choose!")

    if wavenum_min is None and wavenum_max is None:
        assert wavelength_max is not None
        assert wavelength_min is not None
        wavenum_min = nm2cm(wavelength_max)
        wavenum_max = nm2cm(wavelength_min)
    else:
        assert wavenum_min is not None
        assert wavenum_max is not None

    # ... temperatures

    if Tgas is None and Trot is None:
        raise ValueError(
            "Choose either Tgas (equilibrium) or Tvib / Trot (non equilibrium)"
        )

    if Tvib is None and Trot is not None or Tvib is not None and Trot is None:
        raise ValueError("Choose both Tvib and Trot")

    # ... others
    if databank is None:
        raise ValueError("Give a databank name")

    if not "save_memory" in kwargs:
        # no need to save intermediary results as
        # factory is used once only
        kwargs["save_memory"] = True

    if "chunksize" in kwargs:
        raise DeprecationWarning("use optimization= instead of chunksize=")

    def _is_at_equilibrium():
        try:
            assert Tvib is None or Tvib == Tgas
            assert Trot is None or Trot == Tgas
            assert overpopulation is None
            if "self_absorption" in kwargs:
                assert kwargs["self_absorption"]  # == True
            return True
        except AssertionError:
            return False

    _equilibrium = _is_at_equilibrium()

    # which columns to keep when loading line database
    if kwargs["save_memory"] >= 2 and _equilibrium:
        drop_columns = "all"
    else:
        drop_columns = "auto"

    # Run calculations
    sf = SpectrumFactory(
        wavenum_min,
        wavenum_max,
        medium=medium,
        molecule=molecule,
        isotope=isotope,
        pressure=pressure,
        wstep=wstep,
        broadening_max_width=broadening_max_width,
        db_use_cached=use_cached,
        verbose=verbose,
        optimization=optimization,
        **kwargs
    )
    if databank == "fetch":  # mode to get databank without relying on  Line databases
        if _equilibrium:
            # Get line database from HITRAN
            # and partition functions from HAPI
            sf.fetch_databank(
                source="astroquery",
                format="hitran",
                parfuncfmt="hapi",  # use HAPI partition functions for equilibrium
                levelsfmt=None,  # no need to load energies
            )
        else:
            # Also get line database from HITRAN, and calculate partition functions
            # with energy levels from built-in constants (not all molecules
            # are supported!)
            sf.fetch_databank(
                source="astroquery",
                format="hitran",
                parfuncfmt="hapi",  # use HAPI partition functions for equilibrium
                levelsfmt="radis",  # built-in spectroscopic constants
            )
    elif exists(databank):
        # Guess format
        if databank.endswith(".par"):
            if verbose:
                print("Infered {0} is a HITRAN-format file.".format(databank))
            # If non-equilibrium we'll also need to load the energy levels.
            if _equilibrium:
                # Get partition functions from HAPI
                sf.load_databank(
                    path=databank,
                    format="hitran",
                    parfuncfmt="hapi",  # use HAPI partition functions for equilibrium
                    levelsfmt=None,  # no need to load energies
                    drop_columns=drop_columns,
                )
            else:
                # calculate partition functions with energy levels from built-in
                # constants (not all molecules are supported!)
                sf.load_databank(
                    path=databank,
                    format="hitran",
                    parfuncfmt="hapi",  # use HAPI partition functions for equilibrium
                    levelsfmt="radis",  # built-in spectroscopic constants
                    drop_columns=drop_columns,
                )
        elif databank.endswith(".npy"):
            if verbose:
                print("Infered {0} is a NPY-format file".format(databank))

            if _equilibrium:
                sf.load_databank(
                    path=databank,
                    format="cdsd-hitemp",
                    parfuncfmt="hapi",
                    levelsfmt=None,
                    buffer="npy",
                )
            else:
                raise (
                    AttributeError(
                        "Non equilibirum spectra calculation not yet supported with npy databank"
                    )
                )
        else:
            raise ValueError(
                "Couldnt infer the format of the line database file: {0}. ".format(
                    databank
                )
                + "Create a user-defined database in your ~/.radis file "
                + "and define the format there. More information on "
                + "https://radis.readthedocs.io/en/latest/lbl/lbl.html#configuration-file"
            )

    else:  # manual mode: get from user-defined line databases defined in ~/.radis
        sf.load_databank(
            databank,
            load_energies=not _equilibrium,  # no need to load/calculate energies at eq.
            drop_columns=drop_columns,
        )

    #    # Get optimisation strategies
    #    if lineshape_optimization == 'auto':        # NotImplemented: finally we use DLM all the time as default.
    #        if len(sf.df0) > 1e5:
    #            lineshape_optimization = 'DLM'
    #        else:
    #            lineshape_optimization = None
    #        sf.params['chunksize'] = lineshape_optimization

    # Use the standard eq_spectrum / non_eq_spectrum functions
    if _equilibrium:
        if mode == "cpu":
            s = sf.eq_spectrum(
                Tgas=Tgas,
                mole_fraction=mole_fraction,
                path_length=path_length,
                name=name,
            )
        else:
            s = sf.eq_spectrum_gpu(
                Tgas=Tgas,
                mole_fraction=mole_fraction,
                pressure=pressure,
                path_length=path_length,
                name=name,
            )
    else:
        s = sf.non_eq_spectrum(
            Tvib=Tvib,
            Trot=Trot,
            Ttrans=Tgas,
            overpopulation=overpopulation,
            mole_fraction=mole_fraction,
            path_length=path_length,
            name=name,
        )

    return s


# --------------------------
if __name__ == "__main__":

    from radis.test.lbl.test_calc import _run_testcases

    print(_run_testcases(verbose=True))
