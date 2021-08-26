# -*- coding: utf-8 -*-
"""
Define warnings for radiation code, and how to deal with them

----------


Main warning classes :
- :py:class:`~radis.misc.warning.AccuracyWarning`
- :py:class:`~radis.misc.warning.PerformanceWarning`
- default ``UserWarning``


"""


import warnings

from radis.misc.printer import printr

# %% Spectrum warnings / errors
# -----------------------------


class SlitDispersionWarning(UserWarning):
    """Warning trigger if Slit dispersion is too large."""


# %% Spectrum Factory warnings / errors
# -------------------------------------


# Main warning classes :
# - AccuracyWarning
# - PerformanceWarning
# - default UserWarning


class AccuracyError(ValueError):
    """Output spectrum is not valid"""

    pass


class AccuracyWarning(UserWarning):
    """Warning triggered when it seems accuracy is low."""

    pass


class PerformanceWarning(UserWarning):
    """Warning triggered when it seems computation parameters are not
    optimized."""

    pass


# Other warnings


class OutOfBoundError(ValueError):
    pass


class OutOfBoundWarning(UserWarning):
    """Out of bound (for partition functions)"""

    pass


# Errors


class EmptyDatabaseError(ValueError):
    pass


# Warnings


class GaussianBroadeningWarning(AccuracyWarning):
    pass


class CollisionalBroadeningWarning(AccuracyWarning):
    pass


class VoigtBroadeningWarning(AccuracyWarning):
    pass


class MemoryUsageWarning(PerformanceWarning):
    pass


class EmptyDatabaseWarning(UserWarning):
    """Trigger a warning if Line database is empty in the range considered."""

    pass


class OutOfRangeLinesWarning(UserWarning):
    """Trigger a warning if out of range neighbouring lines, that could have an
    effect on the spectrume due to their broadening, cannot be found in the
    database."""

    pass


class HighTemperatureWarning(UserWarning):
    """Warning triggered when the Line database seems inappropriate for the
    temperatures considered."""

    pass


class NegativeEnergiesWarning(UserWarning):
    pass


class MissingSelfBroadeningTdepWarning(UserWarning):
    """Self broadening temperature dependance coefficient is missing in Line Database.

    Usually, use Air broadening temperature dependance coefficient instead. See
    :py:meth:`~radis.lbl.broadening.BroadenFactory._add_collisional_broadening_HWHM`
    """

    pass


class MissingSelfBroadeningWarning(UserWarning):
    """Self broadening tabulated width is missing in Line Database.

    Usually, use Air broadening tabulated width instead. See
    :py:meth:`~radis.lbl.broadening.BroadenFactory._add_collisional_broadening_HWHM`
    """

    pass


class MissingPressureShiftWarning(UserWarning):
    """Pressure-shift coefficient is missing in Line Database."""

    # TODO : add docstring link to references of line database columns.
    pass


class LinestrengthCutoffWarning(AccuracyWarning):
    """Warning triggered when the cumulated linestrength after intensity cutoff
    has changed too much."""

    pass


class InputConditionsWarning(UserWarning):
    """Warning triggered when Spectrum input conditions are suspicious."""

    pass


class DeprecatedFileWarning(DeprecationWarning):
    """Warning triggered when the cached file was generated in a previous version of radis"""

    pass


class IrrelevantFileWarning(PerformanceWarning):
    """Warning triggered when the cached file is irrelevant for the current calcul"""

    pass


class MissingReferenceWarning(UserWarning):
    """Warning triggered when some algorithm / database is missing the bibliographic
    data used by :py:meth:`~radis.spectrum.spectrum.Spectrum.cite`"""

    pass


# %% Config file warnings & errors


class DatabaseAlreadyExists(KeyError):
    pass


class DatabaseNotFoundError(FileNotFoundError):
    """Warning triggered when path does not exist"""

    pass


# @dev: list all your custom warnings below so they are handled by RADIS user params.
WarningClasses = {
    "default": UserWarning,
    "AccuracyWarning": AccuracyWarning,
    "PerformanceWarning": PerformanceWarning,
    "AccuracyError": AccuracyError,
    "SlitDispersionWarning": SlitDispersionWarning,
    "GaussianBroadeningWarning": GaussianBroadeningWarning,
    "CollisionalBroadeningWarning": CollisionalBroadeningWarning,
    "VoigtBroadeningWarning": VoigtBroadeningWarning,
    "MemoryUsageWarning": MemoryUsageWarning,
    "EmptyDatabaseWarning": EmptyDatabaseWarning,
    "DatabaseAlreadyExists": DatabaseAlreadyExists,
    "DatabaseNotFoundError": DatabaseNotFoundError,
    "OutOfRangeLinesWarning": OutOfRangeLinesWarning,
    "HighTemperatureWarning": HighTemperatureWarning,
    "NegativeEnergiesWarning": NegativeEnergiesWarning,
    "MissingSelfBroadeningTdepWarning": MissingSelfBroadeningTdepWarning,
    "MissingSelfBroadeningWarning": MissingSelfBroadeningWarning,
    "MissingPressureShiftWarning": MissingPressureShiftWarning,
    "LinestrengthCutoffWarning": LinestrengthCutoffWarning,
    "InputConditionsWarning": InputConditionsWarning,
    "DeprecatedFileWarning": DeprecatedFileWarning,
    "IrrelevantFileWarning": IrrelevantFileWarning,
    "OutOfBoundWarning": OutOfBoundWarning,
    "MissingReferenceWarning": MissingReferenceWarning,
}
""" dict: warnings used in RADIS Spectrum calculations.

Setup individual warnings. Value of keys can be:
- 'warning' (default: just trigger a warning)
- 'error' (raises an error on this warning)
- 'ignore'  (do nothing)

The key self.warnings['default'] will set the warning behavior for all
other warnings

You can selectively activate them at runtime by setting the warnings attribute of
:class:`radis.lbl.factory.SpectrumFactory`

See Also
--------

:py:data:`~radis.misc.warning.default_warning_status`
"""
default_warning_status = {
    "default": "warn",  # default
    "AccuracyWarning": "warn",
    "AccuracyError": "error",
    "PerformanceWarning": "warn",
    "SlitDispersionWarning": "warn",
    "GaussianBroadeningWarning": "once",  # once per Spectrum calculation
    "CollisionalBroadeningWarning": "once",  # once per Spectrum calculation
    "VoigtBroadeningWarning": "once",  # once per Spectrum calculation
    # see also self.misc.warning_broadening_threshold for the treshold value
    "LinestrengthCutoffWarning": "warn",
    "MemoryUsageWarning": "warn",
    "EmptyDatabaseWarning": "warn",
    "DatabaseAlreadyExists": "error",
    "DatabaseNotFoundError": "error",
    "OutOfRangeLinesWarning": "warn",
    "HighTemperatureWarning": "warn",
    "NegativeEnergiesWarning": "warn",  # warning if negative energies in database
    # warning if self-broadening abs coefficnet missing (Air is used instead)
    "MissingSelfBroadeningTdepWarning": "warn",
    "MissingSelfBroadeningWarning": "warn",
    "MissingPressureShiftWarning": "warn",
    "InputConditionsWarning": "warn",
    "DeprecatedFileWarning": "warn",
    "IrrelevantFileWarning": "warn",
    "OutOfBoundWarning": "warn",
    "MissingReferenceWarning": "warn",
}
""" dict: default status of warnings used in RADIS Spectrum calculations.

Value of keys can be:

- ``'warning'`` (default: just trigger a warning)
- ``'error'`` (raises an error on this warning)
- ``'ignore'``  (do nothing)

The key self.warnings['default'] will set the warning behavior for all
other warnings. All warnings can be disabled by setting the SpectrumFactory
:py:attr:`~radis.lbl.loader.DatabankLoader.warnings` attribute to ``False``.

See Also
--------

:py:data:`~radis.misc.warning.WarningClasses`,
:py:func:`~radis.misc.warning.reset_warnings`

"""


def reset_warnings(status):
    """Reactivate warnings that are set 'once' per session in the Factory
    (unless all warnings have been set to False)

    Parameters
    ----------
    status: dict
        dictionary of Warnings with associated status
    """

    if status == False:
        return

    for k, v in status.items():
        if v == "once":
            WarningType = WarningClasses[k]
            warnings.simplefilter("default", WarningType)


def warn(message, category="default", status={}):
    """Trigger a warning, an error or just ignore based on the value defined in
    the :py:attr:`~radis.lbl.loader.DatabankLoader.warnings` dictionary.

    The warnings can thus be deactivated selectively by setting the SpectrumFactory
    :attr:`~radis.lbl.loader.DatabankLoader.warnings` attribute. All warnings
    can be disabled by setting it to ``False``.

    Parameters
    ----------
    message: str
        what to print
    category: str
        one of the keys of self.warnings.
    status: dict
        status for all warning categories. Can be one of ``'warn'``, ``'ignore'``,
        ``'print'``, ``'error'``

    Examples
    --------
    ::

        if not ((df.Erotu > tol).all() and (df.Erotl > tol).all()):
            warn(
                "There are negative rotational energies in the database",
                "NegativeEnergiesWarning",
            )

    """
    # TODO (refactor): make it possible to run warn(NegativeEnergiesWarning("message"))
    # instead of warn("message, "NegativeEnergiesWarning")
    # ex :
    # if isinstance(message, Warning):
    #     etc.

    if status == False:
        return

    action = status[category]

    WarningType = WarningClasses[category]

    if action in "warn":
        warnings.warn(WarningType(message))
    elif action == "once":
        warnings.warn(WarningType(message))
        # keep 'once' but ignore WarningType with simplefilters
        warnings.simplefilter("ignore", WarningType)
    elif action == "ignore":
        pass
    elif action == "print":  # just print the message, in red
        printr(message)
    elif action == "error":
        raise WarningType(message)
    else:
        raise ValueError("Unexpected action for warning: {0}".format(action))


#%% Tests (on module load)


# ... test warnings are well defined
for k in default_warning_status.keys():
    assert k in WarningClasses

# ... and reciprocally, but they all have a default value
for k in WarningClasses.keys():
    assert k in default_warning_status
