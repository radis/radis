# -*- coding: utf-8 -*-
"""
Functions and constants used in :class:`~radis.spectrum.spectrum.Spectrum` 
object

-------------------------------------------------------------------------------


"""

from __future__ import print_function, absolute_import, division, unicode_literals
from radis.misc.basics import partition

# %% Definitions

# Waverange, wavespaces
WAVENUM_UNITS = ['cm', 'cm-1', 'cm_1', 'wavenumber']
WAVELEN_UNITS = ['nm', 'wavelength']

# Spectral quantities
CONVOLUTED_QUANTITIES = ['radiance', 'transmittance', 'emissivity']
'''list: name of spectral quantities after convolution with slit function'''
NON_CONVOLUTED_QUANTITIES = ['radiance_noslit', 'transmittance_noslit',
                             'emisscoeff', 'emisscoeff_continuum', 
                             'absorbance', 'abscoeff',
                             'abscoeff_continuum', 'emissivity_noslit']
'''list: name of spectral quantities not convolved with slit function'''

# note: it is hardcoded (and needed) that quantities that are convoluted are
# generated from a non convoluted quantity with the same name + _noslit
for _q in CONVOLUTED_QUANTITIES:
    assert _q+'_noslit' in NON_CONVOLUTED_QUANTITIES

# %% Special names when printing conditions

# These parameters are shown below "Physical Conditions" rather than "Computation
# Parameters
PHYSICAL_PARAMS = ['molecule', 'wavenum_max', 'wavenum_min', 'mole_fraction',
                   'isotope', 'state', 'path_length', 'medium', 'self_absorption',
                   'slit_function_base', 'pressure_mbar',
                   'wavelength_min', 'wavelength_max',
                   'Telec', 'Tvib', 'Trot', 'Tgas', 'vib_distribution', 'rot_distribution',
                   'overpopulation', 'thermal_equilibrium']
'''list: physical conditions under which the Spectrum was calculated/measured. 
When printing an object, these parameters are shown below "Physical Conditions" 
If a parameter is not in this list, it is either in "Computation Parameters" 
(non-physical parameters that can have an influence on the Spectrum, e.g, cutoffs
and thresholds) or in "Informative Params" (descriptive parameters that have absolutely no
impact on the spectrum, e.g, number of lines calculated or calculation time)
'''

# Informative parameters
# ... Parameters that should be saved in the Spectrum objects, but ignored when
# ... comparing two spectra.
# ... should be written here only these parameters that cannot affect the
# ... physical result. In particular, all parameters relative to performance
# ... should be added here.
INFORMATIVE_PARAMS = [
    'db_use_cached',
    'chunksize',
    'calculation_time',
    'lines_calculated',
    'lines_cutoff',
    'lines_in_continuum',
    'parallel',
    'Nprocs',
    'Ngroups',
]
''' list: Informative parameters. Parameters that should be saved in the Spectrum 
objects, but ignored when comparing two spectra. Should be written here only 
these parameters that cannot affect the physical result. In particular, all 
parameters relative to performance should be added here.'''

# %% Util functions


def cast_waveunit(unit, force_match=True):
    ''' Standardize unit formats '''
    if unit in WAVELEN_UNITS:
        return 'nm'
    elif unit in WAVENUM_UNITS:
        return 'cm-1'
    elif force_match:
        raise ValueError('Unknown wavespace unit: {0}. Should be one of {1}'.format(unit,
                                                                                    WAVELEN_UNITS+WAVENUM_UNITS))
    else:
        return unit  # dont convert


def make_up(label):
    ''' Cosmetic changes on label, before plot 


    Parameters    
    ----------

    label: str
    '''

    # Improve units
    label = label.replace('cm_1', 'cm-1')
    label = label.replace('cm-1', 'cm$^{-1}$')
    label = label.replace('m2', 'm$^2$')
    label = label.replace('m3', 'm$^3$')
    label = label.replace('I/I0', 'I/I$_\mathrm{0}$')    # transmittance unit

    # Improve text
#    label = label.replace('_noslit', ' (unconvolved)')
    label = label.replace('_noslit', '')     # make no difference between no_slit and slit while plotting
    return label


def print_conditions(conditions, units,
                     phys_param_list=PHYSICAL_PARAMS, info_param_list=INFORMATIVE_PARAMS):
    ''' Print all Spectrum calculation parameters 

    Parameters
    ----------

    phys_param_list: list
        These parameters are shown below "Physical Conditions" rather than "Computation
        Parameters. See :data:`~radis.spectrum.utils.PHYSICAL_PARAMS` for more 
        information. 

    info_param_list: list
        These parameters are shown below "Information" rather than "Computation 
        Parameters. See :data:`~radis.spectrum.utils.INFORMATIVE_PARAMS` for more 
        information.
        
    See Also
    --------
    
    :data:`~radis.spectrum.utils.PHYSICAL_PARAMS`, :data:`~radis.spectrum.utils.INFORMATIVE_PARAMS`

    '''

    def align(a, space=20):
        ''' fix alignement '''
        return a + ' ' * max(1, (space - len(str(a))))

    def print_param(k):
        ''' Special formatting for nicely printing conditions '''
        v_k = conditions[k]
        # Add extra arguments based on arbitrary conditions
        args = []
        if k in units:
            args.append(units[k])
        # ... fill here for other args

        # Special formatting
        try:
            if k in ['wavenum_max_calc', 'wavenum_min_calc', 'wavelength_max', 'wavelength_min',
                     'wavenum_max', 'wavenum_min']:
                v_k_str = '{0:.4f}'.format(v_k)
            elif k in ['lines_calculated', 'lines_in_continuum']:
                # Add comma separator for thousands
                v_k_str = '{0:,d}'.format(v_k)
            else:
                # Default to printing str
                v_k_str = '{0}'.format(v_k)
        except ValueError:
            # Default to printing str
            v_k_str = '{0}'.format(v_k)
    
        # Crop
        if len(v_k_str) > 102:   # cut if too long
            v_k_str = v_k_str[:100] + '...'
        print(' '*2, align(k), v_k_str, *args)

    phys_param, non_phys_param = partition(lambda x: x in phys_param_list,
                                           conditions)

    info_param, non_phys_param = partition(lambda x: x in info_param_list,
                                           non_phys_param)

    print('Physical Conditions')
    print('-'*40)
    for k in sorted(phys_param):
        print_param(k)

    print('Computation Parameters')
    print('-'*40)
    for k in sorted(non_phys_param):
        print_param(k)

    if len(info_param) > 0:
        print('Information')
        print('-'*40)
        for k in sorted(info_param):
            print_param(k)

    print('-'*40)

    # print gas_inp (information on each gas slab) if exists (specifically
    # for Specair output)
    if 'gas_inp' in conditions:
        try:
            for slab in conditions['gas_inp']:
                print('Slab', slab)
                slab.print_conditions()
        except:
            pass

    return None
