# -*- coding: utf-8 -*-
"""

Predefined molecules

See :py:data:`~radis.db.molecules.Molecules`

-------------------------------------------------------------------------------


"""

from __future__ import division, absolute_import, print_function, unicode_literals
from radis.db.classes import ElectronicState
from radis.db.utils import get_dunham_coefficients, get_herzberg_coefficients


# %% Define some commonly used molecules

from radis.phys.convert import eV2cm

# CO
# ----------

# Define with default diatomic constants
CO_X_iso1 = ElectronicState('CO', isotope=1, state='X', term_symbol='1Σ+',
                            spectroscopic_constants=get_dunham_coefficients('CO', 1, 'X1SIG+'),
                            vmax=17,        # max level for Dunham's expansion
                            vmax_morse=48,
                            Ediss=eV2cm(11.16),
                            )
CO_X_iso2 = ElectronicState('CO', isotope=2, state='X', term_symbol='1Σ+',
                            spectroscopic_constants=get_dunham_coefficients('CO', 2, 'X1SIG+'),
                            vmax=17,        # max level for Dunham's expansion
                            vmax_morse=48,
                            Ediss=eV2cm(11.16),
                            )
CO_X_iso3 = ElectronicState('CO', isotope=3, state='X', term_symbol='1Σ+',
                            spectroscopic_constants=get_dunham_coefficients('CO', 3, 'X1SIG+'),
                            vmax=17,        # max level for Dunham's expansion
                            vmax_morse=48,
                            Ediss=eV2cm(11.16),
                            )

# CO2
# ----------

# Generate the energy

#from radis.db.CO2.energies import E_CO2_626_X1SIGg, Ehaj123_CO2_626_X1SIGg
from radis.levels.vibrating_rotor import EvJ_uncoupled_vibrating_rotor, EvJah_uncoupled_vibrating_rotor

#def build_CO2_X(isotope):
#    coeff_dict =get_herzberg_coefficients('CO2', isotope, 'X1SIGu+')
#    # Update standard vibrating rotor formula with CO2 specific generacies and 
#    # Herzberb coeffs
#    def Erovib_CO2(v1, v2, l2, v3, J, offset=True):
#        ''' Rovibrational energies.
#        See :py:func:`~radis.levels.vibrating_rotor.EvJ_uncoupled_vibrating_rotor` 
#        '''
#        return EvJ_uncoupled_vibrating_rotor(v1, v2, l2, v3, J, coeff_dict=coeff_dict,
#                                             gv1=1, gv2=2, gv3=1, offset=offset)
#    def Ehaj_CO2(v1, v2, l2, v3, J, offset=True):
#        ''' Harmonic and anharmonic parts of rovibrational energies.
#        See :py:func:`~radis.levels.vibrating_rotor.EvJah_uncoupled_vibrating_rotor` 
#        '''
#        return EvJah_uncoupled_vibrating_rotor(v1, v2, l2, v3, J, coeff_dict=coeff_dict,
#                                             offset=offset)
#    return ElectronicState('CO2', isotope=1, state='X', term_symbol='1Σu+',
##                            spectroscopic_constants=get_herzberg_coefficients('CO2', 1, 'X1SIGu+'),
#                            Erovib=Erovib_CO2,
#                            Ehaj=Ehaj_CO2,  # for Treanor
#                            Ediss=44600,
#                            )


CO2_X_626 = ElectronicState('CO2', isotope=1, state='X', term_symbol='1Σu+',
                            spectroscopic_constants=get_herzberg_coefficients('CO2', 1, 'X1SIGu+'),
#                            Erovib=EvJah_uncoupled_vibrating_rotor,
#                            Ehaj=Ehaj_CO2,  # for Treanor
                            Ediss=44600,
                            )
CO2_X_636 = ElectronicState('CO2', isotope=2, state='X', term_symbol='1Σu+',
                            spectroscopic_constants=get_herzberg_coefficients('CO2', 2, 'X1SIGu+'),
#                            Erovib=EvJah_uncoupled_vibrating_rotor,
#                            Ehaj=Ehaj_CO2,  # for Treanor
                            Ediss=44600,
                            )
CO2_X_628 = ElectronicState('CO2', isotope=3, state='X', term_symbol='1Σu+',
                            spectroscopic_constants=get_herzberg_coefficients('CO2', 3, 'X1SIGu+'),
#                            Erovib=EvJah_uncoupled_vibrating_rotor,
#                            Ehaj=Ehaj_CO2,  # for Treanor
                            Ediss=44600,
                            )
CO2_X_627 = ElectronicState('CO2', isotope=4, state='X', term_symbol='1Σu+',
                            spectroscopic_constants=get_herzberg_coefficients('CO2', 4, 'X1SIGu+'),
#                            Erovib=EvJah_uncoupled_vibrating_rotor,
#                            Ehaj=Ehaj_CO2,  # for Treanor
                            Ediss=44600,
                            )
#CO2_X_626 = build_CO2_X(1)
#CO2_X_636 = build_CO2_X(2)
#CO2_X_628 = build_CO2_X(3)
#CO2_X_627 = build_CO2_X(4)

#CO2_X_626 = ElectronicState('CO2', isotope=1, state='X', term_symbol='1Σu+',
#                            spectroscopic_constants=get_herzberg_coefficients('CO2', 1, 'X1SIGu+'),
#                            Erovib=EvJ_uncoupled_vibrating_rotor,
#                            Ehaj=EvJah_uncoupled_vibrating_rotor,  # for Treanor
#                            Ediss=44600,
#                            )
#CO2_X_636 = ElectronicState('CO2', isotope=2, state='X', term_symbol='1Σu+',
#                            spectroscopic_constants=get_herzberg_coefficients('CO2', 2, 'X1SIGu+'),
#                            Erovib=EvJ_uncoupled_vibrating_rotor,
#                            Ehaj=EvJah_uncoupled_vibrating_rotor,  # for Treanor
#                            Ediss=44600,
#                            )
#CO2_X_628 = ElectronicState('CO2', isotope=3, state='X', term_symbol='1Σu+',
#                            spectroscopic_constants=get_herzberg_coefficients('CO2', 3, 'X1SIGu+'),
#                            Erovib=EvJ_uncoupled_vibrating_rotor,
#                            Ehaj=EvJah_uncoupled_vibrating_rotor,  # for Treanor
#                            Ediss=44600,
#                            )
#CO2_X_627 = ElectronicState('CO2', isotope=4, state='X', term_symbol='1Σu+',
#                            spectroscopic_constants=get_herzberg_coefficients('CO2', 4, 'X1SIGu+'),
#                            Erovib=EvJ_uncoupled_vibrating_rotor,
#                            Ehaj=EvJah_uncoupled_vibrating_rotor,  # for Treanor
#                            Ediss=44600,
#                            )


# %% Dictionary of predefined molecules
# -----------

#  Molecule  Isotope  ElecState
#         :      :       :
Molecules = {
    'CO':   {
        1:      {
            'X': CO_X_iso1},
        2:      {
            'X': CO_X_iso2},
        3:      {
            'X': CO_X_iso3},
    },
    'CO2':  {
        1:      {
            'X': CO2_X_626},
        2:      {
            'X': CO2_X_636},
        3:      {
            'X': CO2_X_628},
        4:      {
            'X': CO2_X_627},
    },
}
'''dict: list of Electronic states whose energy levels can be calculated with RADIS
with spectroscopic constants. For references refer to the definition of each
molecule. 

See Also
--------

CO:

- :py:data:`~radis.db.molecules.CO_X_iso1`
- :py:data:`~radis.db.molecules.CO_X_iso2`
- :py:data:`~radis.db.molecules.CO_X_iso3`

CO2:
    
- :py:data:`~radis.db.molecules.CO2_X_626`
- :py:data:`~radis.db.molecules.CO2_X_636`
- :py:data:`~radis.db.molecules.CO2_X_628`
- :py:data:`~radis.db.molecules.CO2_X_627`


'''





def getMolecule(molecule, isotope=None, electronic_state=None, verbose=True):
    ''' Get an ElectronicState object in the RADIS :py:data:`~radis.db.molecules.Molecules`
    list. 

    Parameters
    ----------

    molecule: str
        molecule name

    isotope: int, or ``None``
        isotope number. if None, only one isotope must exist in database. Else, 
        an error is raised

    electronic_state: str
        if None, only one electronic state must exist in database. Else, an error 
        is raised

    verbose: boolean
        if ``True``, print which electronic state we got
        
    Returns
    -------
    
    state: ElectronicState
        an :py:class:`~radis.db.classes.ElectronicState` object. 

    '''

    # molecule
    try:
        mol = Molecules[molecule]
    except KeyError:
        raise KeyError('{0} is not defined in molecules with built-in '.format(molecule) +
                       'spectroscopic constants. Choose one of: {0}'.format(list(Molecules.keys())))

    # isotope
    if isotope is None:
        if len(list(mol.keys())) != 1:
            raise ValueError('Please precise which isotope among: {0}'.format(
                list(mol.keys())))
        isotope = list(mol.keys())[0]
    try:
        iso = mol[isotope]
    except KeyError:
        raise KeyError('Isotope {0} is not defined for molecule {1}. Choose one of: {2}'.format(
            isotope, molecule, list(mol.keys())))

    # electronic state
    if electronic_state is None:
        if len(list(iso.keys())) != 1:
            raise ValueError('Please choose which electronic state among: {0}'.format(
                list(iso.keys())))
        electronic_state = list(iso.keys())[0]
    try:
        state = iso[electronic_state]
    except KeyError:
        raise KeyError('{0} is not defined for molecule {1}(iso={2}). Choose one of: {3}'.format(
            electronic_state, molecule, isotope, list(mol.keys())))

    # print name
    if verbose:
        print(('Found {0} in RADIS database'.format(state.get_fullname())))

    # Return
    return state


# %% Test


if __name__ == '__main__':
    
    from radis.test.db.test_molecules import _run_testcases
    print(('Testing molecules.py', _run_testcases()))
