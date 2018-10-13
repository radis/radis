# -*- coding: utf-8 -*-
"""
Define warnings for radiation code, and how to deal with them

----------


"""

from __future__ import absolute_import, unicode_literals, print_function, division

import warnings
from radis.misc.printer import printr

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
    ''' Trigger a warning if Line database is empty in the range considered '''
    pass


class NegativeEnergiesWarning(UserWarning):
    pass


class MissingSelfBroadeningWarning(UserWarning):
    ''' Self broadening is missing in Line Database. Usually, use Air broadening
    instead '''
    pass


class LinestrengthCutoffWarning(UserWarning):
    ''' Warning triggered when the cumulated linestrength after intensity cutoff
    has changed too much'''
    pass


class InputConditionsWarning(UserWarning):
    ''' Warning triggered when Spectrum input conditions are suspicious '''
    pass

class PerformanceWarning(UserWarning):
    ''' Warning triggered when it seems computation parameters are not optimized'''
    pass


WarningClasses = {'default': UserWarning,
                  'GaussianBroadeningWarning': GaussianBroadeningWarning,
                  'CollisionalBroadeningWarning': CollisionalBroadeningWarning,
                  'VoigtBroadeningWarning': VoigtBroadeningWarning,
                  'MemoryUsageWarning': MemoryUsageWarning,
                  'EmptyDatabaseWarning': EmptyDatabaseWarning,
                  'NegativeEnergiesWarning': NegativeEnergiesWarning,
                  'MissingSelfBroadeningWarning': MissingSelfBroadeningWarning,
                  'LinestrengthCutoffWarning': LinestrengthCutoffWarning,
                  'InputConditionsWarning': InputConditionsWarning,
                  'PerformanceWarning':PerformanceWarning,
                  }
''' dict: warnings used in RADIS Spectrum calculations.

You can selectively activate them by setting the warnings attribute of 
:class:`radis.lbl.factory.SpectrumFactory` 
'''

# Setup individual warnings. Value of keys can be:
# - 'warning' (default: just trigger a warning)
# - 'error' (raises an error on this warning)
# - 'ignore'  (do nothing)
# The key self.warnings['default'] will set the warning behavior for all
# other warnings
default_warning_status = {
    'default': 'warn',          # default
    'GaussianBroadeningWarning': 'once',          # once per Spectrum calculation
    'CollisionalBroadeningWarning': 'once',       # once per Spectrum calculation
    'VoigtBroadeningWarning': 'once',             # once per Spectrum calculation
    # see also self.misc.warning_broadening_threshold for the treshold value
    'LinestrengthCutoffWarning': 'warn',
    'MemoryUsageWarning': 'warn',
    'EmptyDatabaseWarning': 'warn',
    'NegativeEnergiesWarning': 'warn',    # warning if negative energies in database
    # warning if self-broadening abs coefficnet missing (Air is used instead)
    'MissingSelfBroadeningWarning': 'warn',
    'InputConditionsWarning': 'warn',
    'PerformanceWarning': 'warn',
}
''' dict: default status of warnings used in RADIS Spectrum calculations.

Value of keys can be:

- 'warning' (default: just trigger a warning)
- 'error' (raises an error on this warning)
- 'ignore'  (do nothing)

The key self.warnings['default'] will set the warning behavior for all
other warnings
'''


def reset_warnings(status):
    ''' Reactivate warnings that are set 'once' per session in the Factory

    Parameters
    ----------

    status: dict
        dictionary of Warnings with associated status

    '''

    for k, v in status.items():
        if v == 'once':
            WarningType = WarningClasses[k]
            warnings.simplefilter('default', WarningType)


def warn(message, category='default', status={}):
    ''' Trigger a warning, an error or just ignore based on the value defined
    in the :attr:`~radis.lbl.loader.DatabankLoader.warnings` dictionary

    The warnings can thus be deactivated selectively by setting the SpectrumFactory
     :attr:`~radis.lbl.loader.DatabankLoader.warnings` attribute

    Parameters
    ----------

    message: str
        what to print

    category: str
        one of the keys of self.warnings

    status: dict
        status for all warning categories. Can be one of ``'warn'``, ``'ignore'``,
        ``'print'``, ``'error'``

    '''

    action = status[category]

    WarningType = WarningClasses[category]

    if action == 'warn':
        warnings.warn(WarningType(message))
    elif action == 'once':
        warnings.warn(WarningType(message))
        # keep 'once' but ignore WarningType with simplefilters
        warnings.simplefilter('ignore', WarningType)
    elif action == 'ignore':
        pass
    elif action == 'print':  # just print the message, in red
        printr(message)
    elif action == 'error':
        raise WarningType(message)
    else:
        raise ValueError('Unexpected action for warning: {0}'.format(action))


# Tests (on module load)



# ... test warnings are well defined
for k in default_warning_status.keys():
    assert k in WarningClasses
