# -*- coding: utf-8 -*-
"""
Define warnings for radiation code, and how to deal with them

----------


"""

from __future__ import absolute_import, unicode_literals, print_function, division

import warnings
from radis.misc.printer import printr

# %% Spectrum warnings / errors
# -----------------------------


class SlitDispersionWarning(UserWarning):
    """ Warning trigger if Slit dispersion is too large 
    """


# %% Spectrum Factory warnings / errors
# -------------------------------------


class OutOfBoundError(ValueError):
    pass


class OutOfBoundWarning(UserWarning):
    """ Out of bound (for partition functions)"""

    pass


# Errors


class EmptyDatabaseError(ValueError):
    pass


# Warnings


class GaussianBroadeningWarning(UserWarning):
    pass


class CollisionalBroadeningWarning(UserWarning):
    pass


class VoigtBroadeningWarning(UserWarning):
    pass


class MemoryUsageWarning(UserWarning):
    pass


class EmptyDatabaseWarning(UserWarning):
    """ Trigger a warning if Line database is empty in the range considered """

    pass


class OutOfRangeLinesWarning(UserWarning):
    """ Trigger a warning if out of range neighbouring lines, that could have 
    an effect on the spectrume due to their broadening, cannot be found in the 
    database """

    pass


class HighTemperatureWarning(UserWarning):
    """ Warning triggered when the Line database seems inappropriate for the 
    temperatures considered
    """

    pass


class NegativeEnergiesWarning(UserWarning):
    pass


class MissingSelfBroadeningWarning(UserWarning):
    """ Self broadening is missing in Line Database. Usually, use Air broadening
    instead """

    pass


class LinestrengthCutoffWarning(UserWarning):
    """ Warning triggered when the cumulated linestrength after intensity cutoff
    has changed too much"""

    pass


class InputConditionsWarning(UserWarning):
    """ Warning triggered when Spectrum input conditions are suspicious """

    pass


class PerformanceWarning(UserWarning):
    """ Warning triggered when it seems computation parameters are not optimized"""

    pass


WarningClasses = {
    "default": UserWarning,
    "SlitDispersionWarning": SlitDispersionWarning,
    "GaussianBroadeningWarning": GaussianBroadeningWarning,
    "CollisionalBroadeningWarning": CollisionalBroadeningWarning,
    "VoigtBroadeningWarning": VoigtBroadeningWarning,
    "MemoryUsageWarning": MemoryUsageWarning,
    "EmptyDatabaseWarning": EmptyDatabaseWarning,
    "OutOfRangeLinesWarning": OutOfRangeLinesWarning,
    "HighTemperatureWarning": HighTemperatureWarning,
    "NegativeEnergiesWarning": NegativeEnergiesWarning,
    "MissingSelfBroadeningWarning": MissingSelfBroadeningWarning,
    "LinestrengthCutoffWarning": LinestrengthCutoffWarning,
    "InputConditionsWarning": InputConditionsWarning,
    "PerformanceWarning": PerformanceWarning,
    "OutOfBoundWarning": OutOfBoundWarning,
}
""" dict: warnings used in RADIS Spectrum calculations.

You can selectively activate them by setting the warnings attribute of 
:class:`radis.lbl.factory.SpectrumFactory` 

See Also
--------

:py:data:`~radis.misc.warning.default_warning_status` 
"""

# Setup individual warnings. Value of keys can be:
# - 'warning' (default: just trigger a warning)
# - 'error' (raises an error on this warning)
# - 'ignore'  (do nothing)
# The key self.warnings['default'] will set the warning behavior for all
# other warnings
default_warning_status = {
    "default": "warn",  # default
    "SlitDispersionWarning": "warn",
    "GaussianBroadeningWarning": "once",  # once per Spectrum calculation
    "CollisionalBroadeningWarning": "once",  # once per Spectrum calculation
    "VoigtBroadeningWarning": "once",  # once per Spectrum calculation
    # see also self.misc.warning_broadening_threshold for the treshold value
    "LinestrengthCutoffWarning": "warn",
    "MemoryUsageWarning": "warn",
    "EmptyDatabaseWarning": "warn",
    "OutOfRangeLinesWarning": "warn",
    "HighTemperatureWarning": "warn",
    "NegativeEnergiesWarning": "warn",  # warning if negative energies in database
    # warning if self-broadening abs coefficnet missing (Air is used instead)
    "MissingSelfBroadeningWarning": "warn",
    "InputConditionsWarning": "warn",
    "PerformanceWarning": "warn",
    "OutOfBoundWarning": "warn",
}
""" dict: default status of warnings used in RADIS Spectrum calculations.

Value of keys can be:

- 'warning' (default: just trigger a warning)
- 'error' (raises an error on this warning)
- 'ignore'  (do nothing)

The key self.warnings['default'] will set the warning behavior for all
other warnings. All warnings can be disabled by setting the SpectrumFactory
:py:attr:`~radis.lbl.loader.DatabankLoader.warnings` attribute to ``False``.

See Also
--------

:py:data:`~radis.misc.warning.WarningClasses`, 
:py:func:`~radis.misc.warning.reset_warnings`

"""


def reset_warnings(status):
    """ Reactivate warnings that are set 'once' per session in the Factory
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
    """ Trigger a warning, an error or just ignore based on the value defined
    in the :py:attr:`~radis.lbl.loader.DatabankLoader.warnings` dictionary

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

    """

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


# Tests (on module load)


# ... test warnings are well defined
for k in default_warning_status.keys():
    assert k in WarningClasses

# ... and reciprocally, but they all have a default value
for k in WarningClasses.keys():
    assert k in default_warning_status
