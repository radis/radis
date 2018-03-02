# -*- coding: utf-8 -*-
"""
Created on Sat Nov 11 13:30:16 2017

@author: erwan

Functions to update Spectrum with new spectral quantities that can be derived
from existing ones, or rescale path_length or mole_fraction, or add overpopulations

Most of these are binded as methods to the Spectrum class, but stored here to
unload the spectrum.py file 

"""

import numpy as np
from numpy import log as ln
from numpy import inf, exp
from radis.misc.debug import printdbg
from radis.misc.utils import DatabankNotFound
from radis.lbl.equations import calc_radiance
from warnings import warn
from six import string_types

CONVOLUTED_QUANTITIES = ['radiance', 'transmittance', 'emissivity']
NON_CONVOLUTED_QUANTITIES = ['radiance_noslit', 'transmittance_noslit',
                              'emisscoeff', 'absorbance', 'abscoeff',
                              'abscoeff_continuum', 'emissivity_noslit',]

# note: it is hardcoded (and needed) that quantities that are convoluted are 
# generated from a non convoluted quantity with the same name + _noslit
for _q in CONVOLUTED_QUANTITIES:
    assert _q+'_noslit' in NON_CONVOLUTED_QUANTITIES


def _build_update_graph(spec):
    ''' Find dependencies and equivalences between different spectral parameters
    based on the spectrum conditions (equilibrium?, optically thin?, path length given?)
    '''
    
    path_length = 'path_length' in spec.conditions
    slit = ('slit_function' in spec.conditions and 'slit_unit' in spec.conditions
            and 'norm_by' in spec.conditions)
    equilibrium = spec.is_at_equilibrium()
    optically_thin = spec.is_optically_thin()
    
    all_keys = [
             'abscoeff',
             'absorbance',
             'emisscoeff',
             'emissivity_noslit',
             'transmittance',
             'radiance',
             'radiance_noslit',
             'transmittance_noslit',
             ]
    
    # Build edges of relationships
    derivation = {             # {keys, [list of keys]}
        'transmittance_noslit':[['absorbance']],
        'absorbance':[['transmittance_noslit']],
        }
    
    def derives_from(what, *from_keys):
        for k in from_keys:
            if type(k) == str:
                k = [k]
            try:
                derivation[what].append(k)
            except KeyError:
                derivation[what] = [k]
    
    # Build more equivalence relationships if path_length is given
    if path_length:
        derives_from('abscoeff', 'absorbance')
        derives_from('absorbance', 'abscoeff')
        if optically_thin:
            derives_from('radiance_noslit', 'emisscoeff')
            derives_from('emisscoeff', 'radiance_noslit')
        else:
            derives_from('radiance_noslit', ['emisscoeff', 'abscoeff'])
            derives_from('emisscoeff', ['radiance_noslit', 'abscoeff'])

    if slit:
        if __debug__: printdbg('... build_graph: slit given > convoluted keys can be recomputed')
        derives_from('radiance', 'radiance_noslit')
        derives_from('transmittance', 'transmittance_noslit')
        derives_from('emissivity', 'emissivity_noslit')

    if equilibrium:
        if __debug__: printdbg('... build_graph: equilibrium > all keys derive from one')
        # Anything can be recomputed from anything
        for key in all_keys:
            all_but_k = [[k] for k in all_keys if k != key]
            derives_from(key, *all_but_k)
    
    return derivation

def get_redundant(spec):
    ''' Returns a dictionary of all spectral quantities in spectrum and whether
    they are redundant
    
    Use
    -------
    
    redundant = get_redundant(spec)
    
    '''
    
    derivation_graph = _build_update_graph(spec)
    
    ordered_keys = [  # ranked by priority
             'abscoeff',
             'emisscoeff',
             'absorbance',
             'radiance_noslit',
             'transmittance_noslit',
             'emissivity',
             'emissivity_noslit',
             'radiance',
             'transmittance',
             ]
    
    activated = dict().fromkeys(ordered_keys, False)
    for k in spec.get_vars():
        activated[k] = True
    redundant = dict().fromkeys(ordered_keys, False)
    
    # Parse graph
    for key in ordered_keys[::-1]:    # roots
        if key in derivation_graph:
            for from_keys in derivation_graph[key]:
                if all([activated[k] and not redundant[k] for k in from_keys]):
                    redundant[key] = True
                    continue
        else:
            del redundant[key]
            # Should look at dependencies again ??
            
    return redundant

def get_reachable(spec): #, derivation_graph):
    ''' Get the list of all quantities that can be derived from current available
    quantities, based on given spec conditions '''
    
    derivation_graph = _build_update_graph(spec)
    
    ordered_keys = [  # ranked by priority
             'abscoeff',
             'emisscoeff',
             'absorbance',
             'radiance_noslit',
             'transmittance_noslit',
             'emissivity',
             'emissivity_noslit',
             'radiance',
             'transmittance',
             ]
    
#    activated = dict().fromkeys(ordered_keys, False)
    reachable = dict().fromkeys(ordered_keys, False)
    for k in spec.get_vars():
        reachable[k] = True
    
    # Parse graph
    restart = True
    while restart:
        restart = False
        for key in ordered_keys[::-1]:    # roots
            if key in derivation_graph:
                for from_keys in derivation_graph[key]:  # all different ways to compute this value
                    if all([reachable[k] for k in from_keys]):
                        if not reachable[key]:
                            reachable[key] = True
                            restart = True   # status changed -> restart?
        
    return reachable



def update(spec, quantity='all', optically_thin='default', verbose=True):
    ''' Calculate missing quantities that can be derived from the current quantities
    and conditions
    
    e.g: if path_length and emisscoeff are given, radiance_noslit can be recalculated
    if in an optically thin configuration, else if abscoeff is also given 
    
    
    Parameters    
    ----------
    
    spec: Spectrum
    
    quantity: str
        name of the spectral quantity to recompute. If 'same', only the quantities
        in the Spectrum are recomputed. If 'all', then all quantities that can
        be derived are recomputed. Default 'all'. 

    optically_thin: True, False, or 'default'
        determines whether to calculate radiance with or without self absorption.
        If 'default', the value is determined from the self_absorption key
        in Spectrum.conditions. If not given, False is taken. Default 'default'
        Also updates the self_absorption value in conditions (creates it if 
        doesnt exist)
        
    '''

    # Check inputs
    # ------------

    # Get path length
    if 'path_length' in list(spec.conditions.keys()):
        path_length = spec.conditions['path_length']
        true_path_length = True
    else:
        path_length = 1        # some stuff can still be updated.
        true_path_length = False
        
    # Update optically thin
    if optically_thin not in [True, False, 'default']:
        raise ValueError("optically_thin must be one of True, False, 'default'")
    if optically_thin == 'default':
        if 'self_absorption' in list(spec.conditions.keys()):
            optically_thin = not spec.conditions['self_absorption']
        else:
            optically_thin = False
    old_self_absorption = spec.conditions.get('self_absorption')
    if old_self_absorption != (not optically_thin):
        spec.conditions['self_absorption'] = not optically_thin
        if verbose:
            print('self absorption set to:', spec.conditions['self_absorption'])

    initial = spec.get_vars()

    # This is where everything happens:
    # ---------------------------------
    _recalculate(spec, quantity, path_length, path_length, 1, 1,
                 true_path_length=true_path_length,
                 verbose=verbose)
    
    # Output
    # ------
    
    # Get list of new quantities 
    new_q = [k for k in spec.get_vars() if k not in initial]
    if verbose:
        print('New quantities added: {0}'.format(new_q))
        
    # Final checks
    for k in new_q:
        # make sure units exist
        if not k in spec.units:
            raise ValueError('{0} added but unit is unknown'.format(k))
        
# Rescale functions

# ... absorption coefficient
def rescale_abscoeff(spec, rescaled, initial, old_mole_fraction, new_mole_fraction,
                      old_path_length, waveunit, units, extra, true_path_length):
    '''
    
    Parameters    
    ----------
    
    spec: Spectrum
    
    old_path_length: float
        path length in cm
        
    '''
    
    unit = None

    # First get initial abscoeff
    # ---------------------
    if 'abscoeff' in initial:
        _, abscoeff = spec.get('abscoeff', wunit=waveunit)
    elif 'absorbance' in initial and true_path_length: # aka: true path_lengths given
        if __debug__: printdbg('... rescale: abscoeff k1 = A/L1')
        _, A = spec.get('absorbance', wunit=waveunit)
        abscoeff = A/old_path_length                      # recalculate
        unit = 'cm-1'
    elif 'transmittance_noslit' in initial and true_path_length:
        if __debug__: printdbg('... rescale: abscoeff k1 = -ln(T1)/L1')
        # Get abscoeff from transmittance
        _, T1 = spec.get('transmittance_noslit', wunit=waveunit)
        
        # We'll have a problem if the spectrum is optically thick
        b = (T1== 0)  # no transmittance: optically thick mask
        if b.sum()>0:
            msg = "Transmittance is satured. Can't infer abscoeff. Please give absorbance"
            if 'abscoeff' in extra: # cant calculate this one but let it go
                abscoeff = None
                if __debug__: printdbg(msg)
            else:
                raise ValueError(msg)
            
        # Else, let's calculate it
        abscoeff = -ln(T1)/old_path_length         # recalculate
        unit = 'cm-1'
    elif 'abscoeff' in extra: # cant calculate this one but let it go
        abscoeff = None
    else:
        raise ValueError("Can't rescale abscoeff if transmittance_noslit ({0}) ".format(
                            'transmittance_noslit' in initial)+\
                         'or absorbance ({0}), '.format('absorbance' in initial)+\
                         'and true_path_length ({0}) '.format(true_path_length)+\
                         'are not given. Use optically_thin?')
        
    # Then export rescaled value
    # --------------------
    if abscoeff is not None:
        if __debug__: printdbg('... rescale: abscoeff k2 = k1 * N2/N1')
        abscoeff *= new_mole_fraction / old_mole_fraction     # rescale
        rescaled['abscoeff'] = abscoeff
    if unit is not None:
        units['abscoeff'] = unit

    return rescaled, units

# ... all, if equilibrium and abscoeff was rescaled
def _recompute_all_at_equilibrium(spec, rescaled, wavenumber, Tgas, 
                                  new_path_length, true_path_length,
                                  units):
    ''' 
    
    Parameters    
    ----------
    
    rescaled: dict
        abscoeff must be rescaled already 
    '''
    
    def get_unit_radiance():
        return spec.units.get('radiance_noslit', 'mW/cm2/sr/nm')
    
    def get_unit_emisscoeff(unit_radiance):
        if '/cm2' in unit_radiance:
            return unit_radiance.replace('/cm2', '/cm3')
        else:
            return unit_radiance + '/cm'  # will be simplified by Pint afterwards

    assert true_path_length
    
    abscoeff = rescaled['abscoeff']
    path_length = new_path_length
    
    absorbance = abscoeff*path_length

    # Generate output quantities
    transmittance_noslit = exp(-absorbance)
    emissivity_noslit = 1 - transmittance_noslit
    radiance_noslit = calc_radiance(wavenumber, emissivity_noslit, Tgas, 
                                    unit=get_unit_radiance())
    b = (abscoeff==0)  # optically thin mask
    emisscoeff = np.empty_like(abscoeff)
    emisscoeff[b] = radiance_noslit[b]/path_length              # recalculate (opt thin)
    emisscoeff[~b] = radiance_noslit[~b]/(1-transmittance_noslit[~b])*abscoeff[~b]    # recalculate (non opt thin)
    
    # ----------------------------------------------------------------------

    rescaled['absorbance'] = absorbance
    rescaled['transmittance_noslit'] = transmittance_noslit
    rescaled['emissivity_noslit'] = emissivity_noslit
    rescaled['radiance_noslit'] = radiance_noslit
    rescaled['emisscoeff'] = emisscoeff
    
    
    units['abscoeff'] = 'cm_1'
    units['absorbance'] = '-ln(I/I0)'
    units['transmittance_noslit'] = 'I/I0'
    units['emissivity_noslit'] = 'eps'
    units['radiance_noslit'] = get_unit_radiance()
    units['emisscoeff'] = get_unit_emisscoeff(units['radiance_noslit'])
    
    return rescaled, units

# ... emission coefficient
def rescale_emisscoeff(spec, rescaled, initial, old_mole_fraction, new_mole_fraction,
                        old_path_length, optically_thin, waveunit, units, 
                        extra, true_path_length):
    '''
    
    Parameters    
    ----------
    
    spec: Spectrum
    '''
    
    unit = None
    def get_unit(unit_radiance):
        if '/cm2' in unit_radiance:
            return unit_radiance.replace('/cm2', '/cm3')
        else:
            return unit_radiance + '/cm'  # will be simplified by Pint afterwards

    # Firt get initial emisscoeff j1
    # -------------------

    if 'emisscoeff' in initial:
        if __debug__: printdbg('... rescale: emisscoeff j1 = j1')
        _, emisscoeff = spec.get('emisscoeff', wunit=waveunit, Iunit=units['emisscoeff'])
        
    elif 'radiance_noslit' in initial and true_path_length and optically_thin:
        if __debug__: printdbg('... rescale: emisscoeff j1 = I1/L1')
        _, I = spec.get('radiance_noslit', wunit=waveunit, Iunit=units['radiance_noslit'])
        emisscoeff = I/old_path_length   # recalculate
        unit = get_unit(units['radiance_noslit'])
        
    elif ('radiance_noslit' in initial and true_path_length and
          'abscoeff' in initial):
        if __debug__: printdbg('... rescale: emisscoeff j1 = k1*I1/(1-exp(-k1*L1))')
        # get emisscoeff from (initial) abscoeff and (initial) radiance
        _, I = spec.get('radiance_noslit', wunit=waveunit, Iunit=units['radiance_noslit'])
        _, k = spec.get('abscoeff', wunit=waveunit, Iunit=units['abscoeff'])
        
        # Recalculate in the optically thin range (T=1) and elsewhere
        b = (k==0)  # optically thin mask
        emisscoeff = np.empty_like(k)
        # ... optically thin case
        emisscoeff[b] = I[b]/old_path_length              # recalculate (opt thin)
        
        # ... non optically thin case:
        # ... recalculate transmittance from abscoeff
        T_b = exp(-k[~b]*old_path_length) 
        # ... and solve the RTE on an homogeneous slab
        emisscoeff[~b] = k[~b]*I[~b]/(1-T_b)              # recalculate (non opt thin)
        unit = get_unit(units['radiance_noslit'])
        
    elif ('radiance_noslit' in initial and true_path_length and
          'transmittance_noslit' in initial):
        if __debug__: printdbg('... rescale: emisscoeff j1 = k1*I1/(1-T1)')
        # get emisscoeff from (initial) transmittance and (initial) radiance
        _, I = spec.get('radiance_noslit', wunit=waveunit, Iunit=units['radiance_noslit'])
        _, T = spec.get('transmittance_noslit', wunit=waveunit, Iunit=units['transmittance_noslit'])
        
        # Recalculate in the optically thin range (T=1) and elsewhere
        b = (T==1)  # optically thin mask
        emisscoeff = np.empty_like(T)
        # ... optically thin case
        emisscoeff[b] = I[b]/old_path_length              # recalculate (opt thin)
        
        # ... non optically thin case:
        # ... recalculate abscoeff from transmittance
        T_b = T[~b]
        k_b = -ln(T_b)/old_path_length       
        # ... and solve the RTE on an homogeneous slab
        emisscoeff[~b] = k_b*I[~b]/(1-T_b)              # recalculate (non opt thin)
        unit = get_unit(units['radiance_noslit'])
        
    else:
        if optically_thin:
            msg = "Can't calculate emisscoeff if true path_length ({0})".format(true_path_length)+\
                  "and initial radiance_noslit ({0}) are not all given".format(
                          'radiance_noslit' in initial)
            if 'emisscoeff' in extra: # cant calculate this one but let it go
                emisscoeff = None
                if __debug__: printdbg(msg)
            else:
                raise ValueError(msg)
        else:
            msg = "Trying to get the emission coefficient (emisscoeff) in non optically "+\
                  "thin case. True path_length ({0}), radiance_noslit ({1}) ".format(
                          true_path_length, 'radiance_noslit' in initial)+\
                  "and abscoeff ({0}) are needed but not all given. ".format(
                          'abscoeff' in initial)+\
                  "Try optically_thin? See known Spectrum conditions with "+\
                  "print(Spectrum)"
            if 'emisscoeff' in extra: # cant calculate this one but let it go
                emisscoeff = None
                if __debug__: printdbg(msg)
            else:
                raise ValueError(msg)
            
    # Then rescale and export
    # -----------
    if emisscoeff is not None:
        if __debug__: printdbg('... rescale: emisscoeff j2 = j1 * N2/N1')
        # Now rescale for mole fractions
        emisscoeff *= new_mole_fraction / old_mole_fraction   # rescale
        rescaled['emisscoeff'] = emisscoeff
    if unit is not None:
        units['emisscoeff'] = unit

    return rescaled, units

# ... absorbance
def rescale_absorbance(spec, rescaled, initial, old_mole_fraction, new_mole_fraction,
                        old_path_length, new_path_length, waveunit, units, 
                        extra, true_path_length):
    '''
    
    Parameters    
    ----------
    
    spec: Spectrum
    '''
    
    unit = None

    if 'absorbance' in initial:
        if __debug__: printdbg('... rescale: absorbance A2 = A1*N2/N1*L2/L1')
        _, absorbance = spec.get('absorbance', wunit=waveunit, Iunit=units['absorbance'])
        absorbance *= new_mole_fraction / old_mole_fraction   # rescale
        absorbance *= new_path_length / old_path_length       # rescale
    elif 'abscoeff' in rescaled and true_path_length:
        if __debug__: printdbg('... rescale: absorbance A2 = j2*L2')
        abscoeff = rescaled['abscoeff']        # mole_fraction already scaled
        absorbance = abscoeff*new_path_length               # calculate
        unit = '-ln(I/I0)'
    else:
        msg = 'Cant recalculate absorbance if absoeff and true path_length are not given'
        if 'absorbance' in extra:  # cant calculate this one but let it go
            absorbance = None
            if __debug__: printdbg(msg)
        else:
            raise ValueError(msg)
            
    # Export rescaled value
    if absorbance is not None:
        rescaled['absorbance'] = absorbance
    if unit is not None:
        units['absorbance'] = unit

    return rescaled, units

# ... transmittance
def rescale_transmittance_noslit(spec, rescaled, initial, old_mole_fraction, new_mole_fraction,
                                  old_path_length, new_path_length, waveunit, units, extra):
    '''
    
    Parameters    
    ----------
    
    spec: Spectrum
    '''
    
    unit = None    
    def get_unit():
        return 'I/I0'

    if 'absorbance' in rescaled:
        if __debug__: printdbg('... rescale: transmittance_noslit T2 = exp(-A2)')
        absorbance = rescaled['absorbance'] # N and L already scaled
        transmittance_noslit = exp(-absorbance)           # recalculate
        unit = get_unit()
    elif 'transmittance_noslit' in initial:
        if __debug__: printdbg('... rescale: transmittance_noslit T2 = '+\
                        'exp( ln(T1) * N2/N1 * L2/L1)')
        # get transmittance from initial transmittance
        _, T1 = spec.get('transmittance_noslit', wunit=waveunit, Iunit=units['transmittance_noslit'])
        
        
        # We'll have a problem if the spectrum is optically thick
        b = (T1 == 0)  # optically thick mask
        if b.sum()>0 and (new_mole_fraction < old_mole_fraction or 
                new_path_length < old_path_length):
            # decreasing mole fractions/ path length could increase the transmittance
            # but this information was lost in the saturation
            msg = 'Transmittance is satured. Cant rescale. Please give absorbance'
            if 'transmittance_noslit' in extra: # cant calculate this one but let it go
                transmittance_noslit = None
                if __debug__: printdbg(msg)
            else:
                raise ValueError(msg)
        # Else, just get absorbance
        absorbance = -ln(T1)
        absorbance *= new_mole_fraction / old_mole_fraction     # rescale
        absorbance *= new_path_length / old_path_length         # rescale
        transmittance_noslit = exp(-absorbance)
    else:
        msg = 'Missing data to rescale transmittance'
        if 'transmittance_noslit' in extra: # cant calculate this one but let it go
            transmittance_noslit = None
            if __debug__: printdbg(msg)
        else:
            raise ValueError(msg)
            
    # Export rescaled value
    if transmittance_noslit is not None:
        rescaled['transmittance_noslit'] = transmittance_noslit
    if unit is not None:
        units['transmittance_noslit'] = unit

    return rescaled, units

# ... radiance_noslit
def rescale_radiance_noslit(spec, rescaled, initial, old_mole_fraction, new_mole_fraction,
                             old_path_length, new_path_length, optically_thin, waveunit,
                             units, extra, true_path_length):
    '''
    
    Parameters    
    ----------
    
    spec: Spectrum
    '''
    
    unit = None
    def get_unit(unit_emisscoeff):
        if '/cm3' in unit_emisscoeff:
            return unit_emisscoeff.replace('/cm3', '/cm2')
        else:
            return unit_emisscoeff + '*cm'

    if 'emisscoeff' in rescaled and true_path_length and optically_thin:
        if __debug__: printdbg('... rescale: radiance_noslit I2 = j2 * L2 '+\
                        '(optically thin)')
        emisscoeff = rescaled['emisscoeff']    # mole_fraction already scaled
        radiance_noslit = emisscoeff*new_path_length      # recalculate
        unit = get_unit(units['emisscoeff'])
    elif ('emisscoeff' in rescaled and 'transmittance_noslit' in rescaled
          and 'abscoeff' in rescaled and true_path_length and not optically_thin): # not optically thin
        if __debug__: printdbg('... rescale: radiance_noslit I2 = j2*(1-T2)/k2')
        emisscoeff = rescaled['emisscoeff']    # mole_fraction already scaled
        abscoeff = rescaled['abscoeff']        # mole_fraction already scaled
        transmittance_noslit = rescaled['transmittance_noslit']  # mole_fraction, path_length already scaled
        b = (abscoeff == 0)  # optically thin mask
        radiance_noslit = np.zeros_like(emisscoeff)         # calculate
        radiance_noslit[~b] = emisscoeff[~b]/abscoeff[~b]*(1-transmittance_noslit[~b])
        radiance_noslit[b] = emisscoeff[b]*new_path_length # optically thin limit
        unit = get_unit(units['emisscoeff'])
    elif 'radiance_noslit' in initial and optically_thin:
        if __debug__: printdbg('... rescale: radiance_noslit I2 = I1*N2/N1*L2/L1 '+\
                        '(optically thin)')
        _, radiance_noslit = spec.get('radiance_noslit', wunit=waveunit, Iunit=units['radiance_noslit'])
        radiance_noslit *= new_mole_fraction / old_mole_fraction    # rescale
        radiance_noslit *= new_path_length / old_path_length        # rescale
    else:
        if optically_thin:
            msg = 'Missing data to rescale radiance_noslit in '+\
                             'optically thin mode. You need at least initial '+\
                             'radiance_noslit ({0}), or scaled emission coefficient ({1}) '.format(
                                     'radiance_noslit' in initial, 'emisscoeff' in rescaled)+\
                             'and true path length ({0}).'.format(true_path_length)
            if 'radiance_noslit' in extra: # cant calculate this one but let it go
                radiance_noslit = None
                if __debug__: printdbg(msg)
            else:
                raise ValueError(msg)
        else:
            msg = 'Missing data to recalculate radiance_noslit. '+\
                             'Try in optically thin mode'
            if 'radiance_noslit' in extra:
                radiance_noslit = None
                if __debug__: printdbg(msg)
            else:
                raise ValueError(msg)
            
    # Export rescaled value
    if radiance_noslit is not None:
        rescaled['radiance_noslit'] = radiance_noslit
    if unit is not None:
        units['radiance_noslit'] = unit

    return rescaled, units

# ... emissivity_noslit
def rescale_emissivity_noslit(spec, rescaled, units, extra, true_path_length):
    '''
    
    Parameters    
    ----------
    
    spec: Spectrum
    '''

    if 'transmittance_noslit' in rescaled:
        if __debug__: printdbg('... rescale: emissivity_noslit e2 = 1 - T2')
        T2 = rescaled['transmittance_noslit']    # transmittivity already scaled
        emissivity_noslit = 1 - T2                  # recalculate
    else:
        msg = 'transmittance_noslit needed to recompute emissivity_noslit'
        if 'emissivity_noslit' in extra: # cant calculate this one but let it go
            emissivity_noslit = None
            if __debug__: printdbg(msg)
        else:
            raise ValueError(msg)

    # Export rescaled value
    if emissivity_noslit is not None:
        rescaled['emissivity_noslit'] = emissivity_noslit
        units['emissivity_noslit'] = 'eps'

    return rescaled, units

def _recalculate(spec, quantity, new_path_length, old_path_length,
                 new_mole_fraction, old_mole_fraction,
                 true_path_length=True, verbose=True):
    ''' General function to recalculate missing quantities. Used in rescale
    and update

    
    Parameters    
    ----------

    spec: Spectrum
    
    quantity: str
        name of the spectral quantity to recompute. If 'same', only the quantities
        in the Spectrum are recomputed. If 'all', then all quantities that can
        be derived are recomputed. 

    true_path_length: boolean
        if False, only relative rescaling (new/old) is allowed. For instance,
        when you dont know the true path_lenth, rescaling absorbance
        with *= new_length/old_length is fine, but abscoeff*new_length is not
        Default True
    '''
    
    optically_thin = spec.is_optically_thin()
    if __debug__: printdbg('... rescale: optically_thin: {0}'.format(optically_thin))

    # Check inputs
    assert quantity in CONVOLUTED_QUANTITIES + NON_CONVOLUTED_QUANTITIES + ['all', 'same']

    # Choose which values to recompute (and store them in the list wanted)
    # ----------
    initial = spec.get_vars()               # quantities initialy in spectrum
    if quantity == 'all':                   # quantities to recompute
        wanted = list(initial)
        greedy = True
    elif quantity == 'same':
        wanted = list(initial)
        greedy = False
    elif type(quantity) in string_types:
        wanted = [quantity]
        greedy = False
    else:
        raise ValueError('unexpected type for quantity: expected str, got '+\
                         '{0} ({1})'.format(quantity, type(quantity)))
    rescaled = {}   # quantities rescaled
    
    # list of quantities that are needed to recomputed what we want
    recompute = wanted        
    if 'radiance' in recompute:
        recompute.append('radiance_noslit')   # technically it's possible without radiance_noslit if in optically_thin mode
    if 'transmittance' in recompute:
        recompute.append('transmittance_noslit')  # same comment
    if 'emissivity' in recompute:
        recompute.append('emissivity_noslit')    # same comment
    if ('radiance_noslit' in initial and not optically_thin):
        recompute.append('emisscoeff') # needed
    if ('absorbance' in initial or 'transmittance_noslit' in initial
                         or 'radiance_noslit' in initial and not optically_thin):
        recompute.append('abscoeff')   # needed

    # choose other parameters to calculate along the way because its easy
    extra = []  
    if greedy:
        # ... let's be greedy: recompute all possible quantities. The list of 
        # all spectral quantities is calculated by parsing a tree in get_reachable
        reachable = get_reachable(spec)
        extra = [k for k, v in reachable.items() if v]
        
    wanted = set(wanted)
    recompute = set(recompute+extra)  # remove duplicates

    # Get units 
    units = spec.units.copy()

    # Recompute!
    # ----------
    waveunit = spec.get_waveunit()  # keep all quantities in same waveunit

    if 'abscoeff' in recompute:
        rescaled, units = rescale_abscoeff(spec, rescaled, initial,     # Todo: remove rescaled = ... Dict is mutable no?
                    old_mole_fraction, new_mole_fraction, old_path_length,
                    waveunit, units, extra, true_path_length)
        
    if (spec.is_at_equilibrium() and 'abscoeff' in recompute and 'Tgas' in spec.conditions
        and greedy):
        wavenumber = spec.get_wavenumber('non_convoluted')
        Tgas = spec.conditions['Tgas']
        rescaled, units = _recompute_all_at_equilibrium(spec, rescaled, wavenumber, Tgas, 
                                                 new_path_length, true_path_length,
                                                 units)
        
    else:
    
        if 'emisscoeff' in recompute:
            rescaled, units = rescale_emisscoeff(spec, rescaled, initial,
                        old_mole_fraction, new_mole_fraction, old_path_length,
                        optically_thin, waveunit, units, extra, true_path_length)
    
        if 'absorbance' in recompute:
            rescaled, units = rescale_absorbance(spec, rescaled, initial,
                        old_mole_fraction, new_mole_fraction, old_path_length,
                        new_path_length, waveunit, units, extra, true_path_length)
    
        if 'transmittance_noslit' in recompute:
            rescaled, units = rescale_transmittance_noslit(spec, rescaled, initial,
                        old_mole_fraction, new_mole_fraction, old_path_length,
                        new_path_length, waveunit, units, extra)
    
        if 'radiance_noslit' in recompute:
            rescaled, units = rescale_radiance_noslit(spec, rescaled, initial,
                        old_mole_fraction, new_mole_fraction, old_path_length,
                        new_path_length, optically_thin, waveunit, units, 
                        extra, true_path_length)
    
        if 'emissivity_noslit' in recompute:
            rescaled, units = rescale_emissivity_noslit(spec, rescaled, units, extra, 
                                                  true_path_length)

    # Save (only) the ones that we want, unless we want everything ('greedy')
    for q in rescaled:
        if q in wanted or greedy:
            if q in NON_CONVOLUTED_QUANTITIES:
                spec._q[q] = rescaled[q]
            else:
                spec._q_conv[q] = rescaled[q]
        
    # Update units
    for k, u in units.items():
        spec.units[k] = u

    # Drop convoluted values
    for q in CONVOLUTED_QUANTITIES:
        if q in list(spec._q_conv.keys()):
            del spec._q_conv[q]

    # Reapply slit if possible
    # TODO: replace with directly convolving with slit stored in conditions
    # TODO: first, add an option to give arrays to apply_slit
    if ('slit_function' in spec.conditions and 'slit_unit' in spec.conditions
        and 'norm_by' in spec.conditions):
        slit_function = spec.conditions['slit_function']
        slit_unit = spec.conditions['slit_unit']
        norm_by = spec.conditions['norm_by']
        try:
            shape = spec.conditions['shape']
        except KeyError:
            shape = 'triangular'
        spec.apply_slit(slit_function=slit_function, unit=slit_unit, shape=shape,
                        norm_by=norm_by, verbose=verbose)            

    # check we didnt miss anyone 
    rescaled_list = list(rescaled)
    for q in spec._q_conv:  # see which one were added by apply_slit
        if not q in rescaled_list:
            rescaled_list.append(q)
    for q in wanted:
        if not q in rescaled_list:
            raise AssertionError('{0} could not be rescaled as wanted. '.format(q)+\
                                 'Could rescale: {0}'.format(rescaled_list))
        
def rescale_path_length(spec, new_path_length, old_path_length=None, force=False):
    ''' Rescale spectrum to new path length. Starts from absorption coefficient
    and emission coefficient, and solves the RTE again for the new path length
    Convoluted values (with slit) are dropped in the process.

    
    Parameters    
    ----------

    spec: Spectrum
    
    new_path_length: float
        new path length

    old_path_length: float, or None
        if None, current path length (conditions['path_length']) is used

        
    Other Parameters
    ----------------
    
    force: boolean
        if False, won't allow rescaling to 0 (not to loose information).
        Default False


    Notes
    -----
    
    Implementation:
    
        To deal with all the input cases, we first make a list of what has to
        be recomputed, and what has to be recalculated

    '''

    # Check inputs
    # ----------
    if old_path_length is not None:
        try:
            if spec.conditions['path_length'] != old_path_length:
                warn('path_length ({0}) doesnt match value given in conditions ({1})'.format(
                        old_path_length, spec.conditions['path_length']))
        except KeyError: # path_length not defined
            pass
    else:
        try:
            old_path_length = spec.conditions['path_length']
        except KeyError:
            raise KeyError('path_length has to be defined in conditions (or use'+\
                            ' `from_path_length`)')

    if new_path_length < 0 and not force:
        raise ValueError('path_length cannot be negative')
    if new_path_length == 0 and not force:
        raise ValueError('Rescaling to 0 will loose information. Choose force '\
                         '= True')
    for q in ['transmittance', 'radiance']:
        qns = q+'_noslit'
        qties = spec.get_vars()
        if q in qties and qns not in qties and not force:
            raise KeyError('Cant rescale {0} if {1} not stored'.format(q, qns)+\
                           ' Use force=True to rescale anyway. {0}'.format(q)+\
                           ' will be deleted')

    # Rescale
    _recalculate(spec, 'same', new_path_length, old_path_length, 1, 1)

    # Update conditions
    spec.conditions['path_length'] = new_path_length


def rescale_mole_fraction(spec, new_mole_fraction, old_mole_fraction=None,
            ignore_warnings=False, force=False):
    ''' Update spectrum with new molar fraction
    Convoluted values (with slit) are dropped in the process.

    
    Parameters    
    ----------

    spec: Spectrum
    
    new_mole_fraction: float
        new mole fraction

    old_mole_fraction: float, or None
        if None, current mole fraction (conditions['mole_fraction']) is used

        
    Other Parameters
    ----------------
    
    force: boolean
        if False, won't allow rescaling to 0 (not to loose information).
        Default False


    Notes
    -----
    
    Implementation:
    
        similar to rescale_path_length() but we have to scale abscoeff & emisscoeff
        Note that this is valid only for small changes in mole fractions. Then,
        the change in line broadening becomes significant


    Todo
    -------

    Add warning when too large rescaling
    '''

    # Check inputs
    # ---------
    if old_mole_fraction is not None:
        try:
            if spec.conditions['mole_fraction'] != old_mole_fraction and not ignore_warnings:
                warn('mole_fraction ({0}) doesnt match value given in conditions ({1})'.format(
                        old_mole_fraction, spec.conditions['mole_fraction']))
        except KeyError: # mole fraction not defined
            pass

    else:
        try:
            old_mole_fraction = spec.conditions['mole_fraction']
        except KeyError:
            raise KeyError('mole_fraction has to be defined in conditions (or use'+\
                            ' `from_mole_fraction`)')

    if new_mole_fraction < 0 and not force:
        raise ValueError('mole_fraction cannot be negative')
    if new_mole_fraction == 0 and not force:
        raise ValueError('Rescaling to 0 will loose information. Choose force '\
                         '= True')

    for q in ['transmittance', 'radiance']:
        qns = q+'_noslit'
        qties = spec.get_vars()
        if q in qties and qns not in qties and not force:
            raise KeyError('Cant rescale {0} if {1} not stored.'.format(q, qns)+\
                           '(you need to rescale before applying the slit again) '+\
                           ' Use force=True to rescale anyway, but {0}'.format(q)+\
                           ' will be deleted')

    # Get path length
    if 'path_length' in list(spec.conditions.keys()):
        path_length = spec.conditions['path_length']
        true_path_length = True
    else:
        path_length = 1
        true_path_length = False

    # Rescale
    _recalculate(spec, 'same', path_length, path_length, new_mole_fraction, old_mole_fraction,
                      true_path_length=true_path_length)

    # Update conditions
    spec.conditions['mole_fraction'] = new_mole_fraction



def _test_compression(verbose=True, warnings=True, *args, **kwargs):
    
    from neq.spec import SpectrumFactory
    
    Tgas = 1500
    sf = SpectrumFactory(
                         wavelength_min=4400,
                         wavelength_max=4800,
    #                     mole_fraction=1,
                         path_length=0.1,
                         mole_fraction=0.01,
                         cutoff=1e-25,
                         wstep = 0.005,
                         isotope=[1],
                         db_use_cached=True,
                         self_absorption=True,
                         verbose=False)
    try:
        sf.load_databank('HITRAN-CO')
    except DatabankNotFound:
        if warnings:
            import sys
            print(sys.exc_info())
            print('Testing spectrum.py: Database not defined: HITRAN-CO \n'+\
                           'Ignoring the test')
        return True

    s1 = sf.non_eq_spectrum(Tgas, Tgas, path_length=0.01)
    redundant = get_redundant(s1)
    if verbose: print(redundant)
    
    return redundant == {'emissivity_noslit': True, 'radiance_noslit': True, 
                         'radiance': True, 'emisscoeff': True, 
                         'transmittance_noslit': True, 'absorbance': True, 
                         'transmittance': True, 'abscoeff': False}

if __name__ == '__main__':

    print('Test rescale.py: ', _test_compression())