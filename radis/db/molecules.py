# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 16:15:03 2017

@author: erwan

Summary
-------

Define Molecule, Isotope and ElectronicState classes 

ElectronicState has:

- Molecular parameters (read from HITRAN molparams.txt)
- Vibrational and rotational energies

Calculate Vibrational and Rotational energies for some database molecules

Molecules
-----------

Shortcut to some molecules ground states are defined at the end of the file 
in the dictionary `Molecules`, with the following structure: 

{molecule_name: {isotope_number: {electronic_state_name: ElectronicState object}}}
        
Import with
    
    >>> from radis.db.molecules import Molecules
    >>> CO2_X = Molecules['CO2'][1]['X']            # 1 for first isotope

Or:

    >>> from radis import getMolecule                # directly
    >>> CO2_2 = getMolecule('CO2', 1, 'X')         # better KeyError messages 


-------------------------------------------------------------------------------


"""

from __future__ import division, absolute_import, print_function, unicode_literals
import inspect
import re
from six import string_types
from radis.io.hitran import get_molecule, get_molecule_identifier
from radis.levels.dunham import Gv, Fv, EvJ
from radis.db.utils import get_dunham_coefficients, get_herzberg_coefficients
from radis.db.conventions import get_convention

# %% Molecule class

class Molecule(object):

    def __init__(self, name, verbose=True):    # no optional kwargs here (final stop)
        ''' Define a new molecule

        Parameters
        ----------

        name: str, int
            molecule name, or HITRAN identifier

        Other Parameters
        ----------------
        
        verbose: boolean
            more blabla
        '''
        
        self.verbose = verbose

        # Get name and integer id
        if isinstance(name, string_types):
            self.name = name

            # Get name without parenthesis (without state) for HITRAN identification
            filtername = re.sub("[\(\[].*?[\)\]]", "", name)

            try:
                self.id = get_molecule_identifier(filtername)
            except KeyError:  # Not an HITRAN molecule
                self.id = None
        elif type(name) == int:
            self.id = name
            self.name = get_molecule(name)
        else:
            raise ValueError('Wrong name type:', name)


# %% Isotope class

class Isotope(Molecule):

    def __init__(self, molecule_name, isotope, isotope_name='', abundance=None,
                 **kwargs):
        ''' Create an isotope 

        Parameters
        ----------

        molecule_name: str, or int
            molecule name or HITRAN identifier

        isotope: int
            isotope identifier, sorted by decreasing abundance (cf HITRAN 
            nomenclature)

        isotope_name: str
            (optional) isotope name. Default ''

        abundance: float
            isotopologue abundance. Default ``None``

        '''

        super(Isotope, self).__init__(
            name=molecule_name, **kwargs)  # initialize Molecule

        # Check input
        if not isinstance(isotope, int):
            raise TypeError('Wrong format for isotope {0}: expected int, got {1}'.format(
                isotope, type(isotope)))

        # Store
        self.iso = isotope
        self.isotope_name = isotope_name
        self.abundance = abundance

# %% ElectronicState class


class ElectronicState(Isotope):
    ''' 
    Notes
    -----
    
    To allow for more flexibility, ElectronicState should also inherit from some 
    built-in Symmetry / species class, Sigma / Pi / Delta etc... where the energy 
    is calculated appropriately.  (ex: Rigid Rotor, Symmetric Top etc.)

    '''
    # TODO  (major code rewritting...)


    def __init__(self, molecule_name, isotope, state, term_symbol='',
                 g_e=None, spectroscopic_constants={}, Erovib=None, Ehaj=None,
                 vmax=None, vmax_morse=None, Jmax=None, Ediss=None, **kwargs):
        ''' Define an electronic state of an isotope of a molecule

        Parameters
        ----------

        molecule_name: str, int
            molecule name, or HITRAN identifier

        isotope: int
            isotope number of the molecule (sorted by abundance on Earth, HITRAN
            convention).

        state: str
            Electronic state name

        term_symbol: str
            Term symbol. Default ''

        g_e: int
            Degeneracy. Default ``None``

        spectroscopic_constants: dict
            spectroscopic constants for given (molecule, isotope, electronic state)
            Expected:
                
                Yij: cm-1
                    rovibrational coefficients in Dunham convention
                
            or
            
                wexe, Be, etc. : cm-1
                    rovibrational coefficients in Herzberg convention
                
            or 
                      
                Te: cm-1
                    electronic energy. Default ``None`` if not given

        Erovib: function
            Rovibrational energy of the molecule. Accessible with, typically, Mol.Erovib(v,J)

        Ehaj: function
            To return vibrational energy in harmonic and anharmonic component, 
            plus rotational energy. Used to build Treanor distributions

        vmax, vmax_morse: int, or ``None``
            maximum vibrational number (required for partition function calculation)
            to be calculated with Dunham expansion, and Morse potential
            If None, number will be infered from dissociation energy

        Jmax: int, or ``None``
            maximum rotational number (required for partition function calculation)
            If None, number will be infered from dissociation energy

        Ediss: cm-1
            dissociation energy. Required for partition function calculation 
            if neither vmax nor Jmax are given

        kwargs: **dict
            forwarded to parent class

        '''

        super(ElectronicState, self).__init__(molecule_name=molecule_name,
                                              isotope=isotope, **kwargs)    # initialize Isotope

        # Initialize electronic state specific parameters
        self.state = state
        self.term_symbol = term_symbol
        self.g_e = g_e

        self._parse_constants(spectroscopic_constants)

        self._assign_E(Erovib)
        self.Ehaj = Ehaj
        self.state = state
        self.vmax = vmax
        self.vmax_morse = vmax_morse
        self.Jmax = Jmax
        self.Ediss = Ediss

    def _parse_constants(self, spectroscopic_constants):
        ''' Parse spectroscopic constants 
        
        Parameters
        ----------
        
        spectroscopic_constants: dict
            spectro constants in Herzberg or Dunham convention
            
        Returns
        -------
        
        None: 
            but constants are stored under self.rovib_constants
        
        '''
        
        import re
        # Clean keys
        # In particular, remove trailing '_cm-1' if given in dict or database 
        rovib_constants = {re.sub('_cm-1$', '', k):v 
                           for (k,v) in spectroscopic_constants.items()}
        
        # Get specific keys
        self.Te = rovib_constants.pop('Te', None)  # default None
        self.re = rovib_constants.pop('re', None)
        
        # Store
        self.rovib_constants = rovib_constants

    # Default method to calculate energy

    def _assign_E(self, Erovib):
        ''' Finds appropriate Electrorovibrational energy function for this molecule ,
        based on which convention is used for spectroscopic coefficients 
        (Dunham, or Herzberg)
        
        Default methods only work for diatomic molecules. For polyatomic molecules,
        you should give the energy function directly (see how CO2 is defined in 
        radis.db.molecules.py) 
        
        Returns
        -------
        
        None:
            but assigns the ``self.Erovib`` function to the molecule
            
        '''
        
        if Erovib is not None:
            # overwrite energy calculation
            self.Erovib = Erovib
            if self.verbose>=2:
                print('{0}: overwritting Energy calculation with {1}'.format(
                        self.get_fullname(), Erovib))
        else:
            c = self.rovib_constants
            if len(c) == 0:
                self.Erovib = self._E_None
                if self.verbose>=2:
                    print('{0}: No rovibconstants found'.format(
                            self.get_fullname()))
            elif get_convention(c) == 'herzberg':
                self.Erovib = self._E_Herzberg
                if self.verbose>=2:
                    print('{0}: using Herzberg coefficients for Energy calculation'.format(
                            self.get_fullname()))
            else:
                self.Erovib = self._E_Dunham
                if self.verbose>=2:
                    print('{0}: using Dunham coefficients for Energy calculation'.format(
                            self.get_fullname()))
                # TODO: reformat these methods as external functions, but keep
                # help

    def _E_None(self, *args, **kwargs):
        ''' Rovibrational Energy '''
        raise NotImplementedError('No spectroscopic coefficients given for '+\
                                  'molecule {0}'.format(self.get_fullname()))
        
    def _E_Dunham(self, v, J, offset=True):
        
        c = self.rovib_constants
        
        if offset:
            ZPE = EvJ(0, 0, **c)
        else:
            ZPE = 0

        return EvJ(v, J, **c) - ZPE
        
    def _E_Herzberg(self, v, J, offset=True):
        ''' Calculate rovibrational energy of molecule 
        
        .. math::
            
            E(v,J) = G(v) + F(v,J)
            
        Works with Herzberg convention for the spectroscopic coefficients used 
        in :func:`~radis.levels.dunham.Gv`, :func:`~radis.levels.dunham.Fv` ::
            
            
        
        Only works for diatomic molecules. Method should be overwritten on molecule
        creation for other molecules 

        Parameters
        ----------

        v: int
            vibrational state 

        J: int
            rotational state

        offset: boolean
            if ``True``, removes energy of v=0,J=0 vibrational level (zero-point-energy)
            Default ``True``

        Returns
        -------

        energy of state in cm-1 
        
        See Also
        --------

        :func:`~radis.levels.dunham.Gv`, :func:`~radis.levels.dunham.Fv`
        
        '''
        

        try:

            # Get from constants (hardcoded, see reference above)
            c = self.rovib_constants
            
            # Vibrational constants
            we = c['we']                 # mandatory
            wexe = c.get('wexe', 0)      # optional: 0 is default value
            weye = c.get('weye', 0)
            weze = c.get('weze', 0)
            weae = c.get('weae', 0)
            webe = c.get('webe', 0)

            # Rotational constants
            Be = c['Be']
            De = c.get('De', 0)

            alpha_e = c.get('alpha_e', 0)
            beta_e = c.get('beta_e', 0)
            gamma_e = c.get('gamma_e', 0)
            delta_e = c.get('delta_e', 0)
            pi_e = c.get('pi_e', 0)

            He = c.get('He', 0)
            eta_e = c.get('eta_e', 0)

            # Energies
#            Te = c['Te']      # electronic energy
            G = Gv(v, we, wexe, weye, weze, weae, webe)   # vibrational energy
            F = Fv(v, J, Be, De, alpha_e, beta_e, gamma_e, delta_e, pi_e=pi_e,
                   He=He, eta_e=eta_e)    # rotational energy

            if offset:
                ZPE = Gv(0, we, wexe, weye, weze)
            else:
                ZPE = 0

#            E = Te + G + F - ZPE    # cm-1   # rovibrational energy only for the moment # TODO
            E = G + F - ZPE    # cm-1

        except KeyError as err:
            raise KeyError('Mandatory spectroscopic constant `{0}` '.format(err.args[0]) +
                           'not defined for electronic state {0}'.format(self.get_fullname()) +
                           '. Check your ElectronicState definition')

        return E

    # %% Informative methods

    def get_Morse_inc(self):
        ''' Get Morse potential energy increment correction for given molecule

        Examples
        --------

        >>> inc = ElecState.get_Morse_inc()
        >>> Delta_E(vi, vi+1) = Delta_E0 - (vi+1-vi0)*inc

        for the energy gap between vibrational level vi and vi+1, 
        where vi0 is the last vibrational level calculated with Dunham 
        expansion, and  and Delta_E0 = DeltaE(vi0-1, vi0)
        
        See Also
        --------
        
        :py:func:`~radis.phys.morse.morse_increment`

        '''
        
        from radis.phys.morse import morse_increment
        
        c = self.rovib_constants

        if get_convention(c) == 'dunham':
        
            # convert to get wexe, weae etc in Herzberg convention
            from radis.db.conventions import herzberg2dunham
            def convert(coef, default=0):
                if coef in c:
                    sign, Yij = herzberg2dunham[coef]
                    return sign*Yij
                else:
                    return default    # default
            
            we = convert('we')             # required
            wexe = convert('wexe', 0)  # optional: 0 is default value
            weye = convert('weye', 0)
            weze = convert('weze', 0)
            weae = convert('weae', 0)
            webe = convert('webe', 0)
            wece = convert('wece', 0)
            
            if we is None:
                raise KeyError('Mandatory spectroscopic constant `we` ' +
                               'not defined for electronic state {0}'.format(self.get_fullname()) +
                               '. Check your ElectronicState definition. Given '+\
                               'constants: {0}'.format(list(c.keys())))

        else:
            try:
                we = c['we']             # required
                wexe = c.get('wexe', 0)  # optional: 0 is default value
                weye = c.get('weye', 0)
                weze = c.get('weze', 0)
                weae = c.get('weae', 0)
                webe = c.get('webe', 0)
                wece = c.get('wece', 0)
                    
            except KeyError as err:
                raise KeyError('Mandatory spectroscopic constant `{0}` '.format(err.args[0]) +
                               'not defined for electronic state {0}'.format(self.get_fullname()) +
                               '. Check your ElectronicState definition. Given '+\
                               'constants: {0}'.format(list(c.keys())))

        return morse_increment(self.Ediss, we, wexe=wexe, weye=weye, weze=weze,
                               weae=weae, webe=webe, wece=wece)

    # Methods to access molecule reference from the console:

    def get_ref(self, var):
        ''' Get references for variable or function `var` 

        Examples
        --------

        >>> from radis.db.molecules import CO
        >>> CO.ref('Erovib')   # get code used to calculate energy E (including  
        >>>                    # corrections to Dunham development constants)

        '''

        att = getattr(self, var)     # raises AttributeError if doesnt exist

        try:
            print(inspect.getsource(att))
        except IOError:
            raise ValueError('Couldnt find energy source code')

    def get_fullname(self):
        return '{0}({1}{2})(iso{3})'.format(self.name, self.state,
                                            self.term_symbol, self.iso)


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

def build_CO2_X(isotope):
    coeff_dict =get_herzberg_coefficients('CO2', isotope, 'X1SIGu+')
    # Update standard vibrating rotor formula with CO2 specific generacies and 
    # Herzberb coeffs
    def Erovib_CO2(v1, v2, l2, v3, J, offset=True):
        ''' Rovibrational energies.
        See :py:func:`~radis.levels.vibrating_rotor.EvJ_uncoupled_vibrating_rotor` 
        '''
        return EvJ_uncoupled_vibrating_rotor(v1, v2, l2, v3, J, coeff_dict=coeff_dict,
                                             gv1=1, gv2=2, gv3=1, offset=offset)
    def Ehaj_CO2(v1, v2, l2, v3, J, offset=True):
        ''' Harmonic and anharmonic parts of rovibrational energies.
        See :py:func:`~radis.levels.vibrating_rotor.EvJah_uncoupled_vibrating_rotor` 
        '''
        return EvJah_uncoupled_vibrating_rotor(v1, v2, l2, v3, J, coeff_dict=coeff_dict,
                                             #gv1=1, gv2=2, gv3=1,   # TODO: look up gv
                                             offset=offset)
    return ElectronicState('CO2', isotope=1, state='X', term_symbol='1Σu+',
#                            spectroscopic_constants=get_herzberg_coefficients('CO2', 1, 'X1SIGu+'),
                            Erovib=Erovib_CO2,
                            Ehaj=Ehaj_CO2,  # for Treanor
                            Ediss=44600,
                            )

CO2_X_626 = build_CO2_X(1)
CO2_X_636 = build_CO2_X(2)
CO2_X_628 = build_CO2_X(3)
CO2_X_627 = build_CO2_X(4)

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

- :class:`~radis.db.molecules.CO_X_iso1`
- :class:`~radis.db.molecules.CO_X_iso2`
- :class:`~radis.db.molecules.CO_X_iso3`

CO2:
    
- :class:`~radis.db.molecules.CO2_X_626`
- :class:`~radis.db.molecules.CO2_X_636`
'''


def getMolecule(molecule, isotope=None, electronic_state=None, verbose=True):
    ''' Get an ElectronicState object the RADIS Molecules list. 

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
        print('Found {0} in RADIS database'.format(state.get_fullname()))

    # Return
    return state


# %% Test


if __name__ == '__main__':
    
    from radis.test.db.test_molecules import _run_testcases
    print('Testing molecules.py', _run_testcases())
