# -*- coding: utf-8 -*-
"""
Summary
-------

Spectrum class holder

Keeps all output data from a SpectrumFactory calculation. Allows to reapply a
different slit function in post-processing, and plot all spectral quantities
with any unit

Routine Listings
----------------

:func:`~radis.spectrum.spectrum.calculated_spectrum`
:func:`~radis.spectrum.spectrum.experimental_spectrum`
:func:`~radis.spectrum.spectrum.transmittance_spectrum`
:func:`~radis.spectrum.spectrum.is_spectrum`


Examples
--------

Typical use::

    from radis calculated_spectrum
    s = calculated_spectrum(w, I, conditions={'case':'previously calculated by ##'})
    s.plot('radiance_noslit')
    s.apply_slit(0.5, shape='triangular')
    s.plot('radiance')

Spectrum objects can be stored, retrieved, rescaled, resamples::

    from radis import load_spec
    s = load_spec('co_calculation.spec')
    s.rescale_path_length(0.5)                  # calculate for new path_length
    s.rescale_mole_fraction(0.02)   # calculate for new mole fraction
    s.resample(w_new)               # resample on new wavespace
    s.store('co_calculation2.spec')

"""

from __future__ import print_function, absolute_import, division, unicode_literals
import matplotlib.pyplot as plt
import numpy as np
from publib import set_style, fix_style
from radis.phys.convert import conv2, cm2nm, nm2cm
from radis.phys.units import Q_, convert_universal
from radis.phys.air import vacuum2air, air2vacuum
from radis.spectrum.utils import (CONVOLUTED_QUANTITIES, NON_CONVOLUTED_QUANTITIES, 
                                  make_up, cast_waveunit, print_conditions)
from radis.spectrum.rescale import update, rescale_path_length, rescale_mole_fraction
#from neq.spec.base import print_conditions
from radis.misc.basics import compare_lists
from radis.misc.arrays import evenly_distributed
from radis.misc.debug import printdbg
from radis.misc.signal import resample
from pint import UndefinedUnitError
from warnings import warn
from numpy import allclose, abs, diff
from copy import deepcopy
from six import string_types

# %% Array-to-Spectrum functions

def calculated_spectrum(w, I, wunit='nm', Iunit='mW/cm2/sr/nm',
                         conditions=None, cond_units=None, populations=None,
                         name=None):
    ''' Convert (w, I) into a Spectrum object that has unit conversion, plotting
    and slit convolution capabilities

    
    Parameters    
    ----------

    w, I: np.array
        wavelength and intensity

    wunit: 'nm', 'cm-1'
        wavespace unit

    Iunit: str
        intensity unit (can be 'counts', 'mW/cm2/sr/nm', etc...). Default
        'mW/cm2/sr/nm' (note that non-convoluted Specair spectra are in 'mW/cm2/sr/µm')


    Other Parameters
    ----------------

    conditions: dict
        (optional) calculation conditions to be stored with Spectrum. Default None

    cond_units: dict
        (optional) calculation conditions units. Default None

    populations: dict
        populations to be stored in Spectrum. Default None

    name: str
        (optional) give a name


    See Also
    --------
    
    :func:`~radis.spectrum.spectrum.transmittance_spectrum`, 
    :func:`~radis.spectrum.spectrum.experimental_spectrum`
    :meth:`~radis.spectrum.spectrum.from_array`
    :meth:`~radis.spectrum.spectrum.from_txt`
    :func:`~radis.tools.database.load_spec`
    
    '''

    return Spectrum.from_array(w, I, 'radiance_noslit', 
                               waveunit=wunit, unit=Iunit,
                               conditions=conditions, cond_units=cond_units,
                               populations=populations,
                               name=name)

def transmittance_spectrum(w, T, wunit='nm', Tunit='I/I0',
                           conditions=None, cond_units=None,
                           name=None):
    ''' Convert (w, I) into a Spectrum object that has unit conversion, plotting
    and slit convolution capabilities

    
    Parameters    
    ----------

    w, I: np.array
        wavelength and transmittance (no slit)

    wunit: 'nm', 'cm-1'
        wavespace unit

    Iunit: str
        intensity unit. Default 'I/I0'


    Other Parameters
    ----------------

    conditions: dict
        (optional) calculation conditions to be stored with Spectrum

    cond_units: dict
        (optional) calculation conditions units

    name: str
        (optional) give a name


    See Also
    --------
    
    :func:`~radis.spectrum.spectrum.calculated_spectrum`, 
    :func:`~radis.spectrum.spectrum.experimental_spectrum`
    :meth:`~radis.spectrum.spectrum.from_array`
    :meth:`~radis.spectrum.spectrum.from_txt`
    :func:`~radis.tools.database.load_spec`
        
    '''

    return Spectrum.from_array(w, T, 'transmittance_noslit', 
                               waveunit=wunit, unit=Tunit,
                               conditions=conditions, cond_units=cond_units,
                               name=name)

def experimental_spectrum(w, I, wunit='nm', Iunit='counts',
                         conditions=None, cond_units=None, name=None):
    ''' Convert (w, I) into a Spectrum object that has unit conversion and plotting
    capabilities. Convolution is not available as the spectrum is assumed to
    be measured experimentally (hence deconvolution of the slit function would
    be required)

    
    Parameters    
    ----------

    w, I: np.array
        wavelength and intensity

    wunit: 'nm', 'cm-1'
        wavespace unit

    Iunit: str
        intensity unit (can be 'counts', 'mW/cm2/sr/nm', etc...). Default
        'counts' (default Winspec output)

    Other Parameters
    ----------------

    conditions: dict
        (optional) calculation conditions to be stored with Spectrum

    cond_units: dict
        (optional) calculation conditions units
        
    name: str
        (optional) give a name


    See Also
    --------
    
    :func:`~radis.spectrum.spectrum.calculated_spectrum`, 
    :func:`~radis.spectrum.spectrum.transmittance_spectrum`, 
    :meth:`~radis.spectrum.spectrum.from_array`
    :meth:`~radis.spectrum.spectrum.from_txt`
    :func:`~radis.tools.database.load_spec`
    
    '''

    return Spectrum.from_array(w, I, 'radiance', 
                               waveunit=wunit, unit=Iunit,
                               conditions=conditions, cond_units=cond_units,
                               name=name)

# %% Spectrum class to hold results )

class Spectrum(object):
    ''' This class holds result from a SpectrumFactory calculation. It can be
    used to plot different quantities a posteriori, or manipulate output units
    (for instance convert a spectral radiance per wavelength units to a
    spectral radiance per wavenumber)

    Parameters
    ----------

    quantities: dict of tuples   {'quantity':(wavenum, quantity)}
        where quantities are spectral quantities (absorbance, radiance, etc.)
        and wavenum is in cm-1
        example:

        >>> {'radiance_noslit:':(wavenum, radiance_noslit),
             'absorbance':(wavenum, absorbance)}

    units: dict
        units for quantities
        
    Other Parameters
    ----------------
    
    Optional parameters:

    conditions: dict
        physical conditions and calculation parameters

    cond_units: dict
        units for conditions

    populations: dict
        a dictionary of all species, and levels. Should be compatible with other
        radiative codes such as Specair output. Suggested format:
        {molecules: {isotopes: {elec state: rovib levels}}}
        e.g::
        
            {'CO2':{1: 'X': df}}   # with df a Pandas Dataframe
        

    lines: pandas Dataframe
        all lines in databank (necessary for LineSurvey). Warning if you want to
        play with the lins content: The signification of columns in `lines` may be
        specific to a database format. Plus, some additional columns may have been
        added by the calculation (e.g: `Ei` and `S` for emission integral and
        linestrength in SpectrumFactory). Refer to the code to know what they mean
        (and their units)

    wavespace: 'nm' or 'cm-1', or None
        define whether wavenumber or wavelength are used in 'quantities' tuples.
        Quantities should be evenly distributed along this space for fast
        convolution with the slit function
        (non-uniform slit function is not implemented anyway... )
        Defaults None (but raises an error if wavespace is not defined in
        conditions neither)


    Other Parameters
    ----------------
    
    name: str, or None
        Give a name to this Spectrum object (helps debugging in multislab
        configurations). Default None

    warnings: boolean
        If True, test if inputs are valid, e.g, spectra are evenly distributed in
        wavelength, and raise a warning if not. Note that this take ~ 3.5 ms for
        a 20k points spectrum, when the rest of the creation process is only 
        ~ 1.8ms (makes it 3 times longer, and can be a problem if hundreds of 
        spectra are created in a row). Default True


    Examples
    --------
    
    Manipulate a Spectrum calculated by RADIS::

        s = calc_spectrum(2125, 2300, Tgas=2000, databank='CDSD')
        s.print_conditions()
        s.plot('absorbance')
        s.line_survey(overlay='absorbance')
        s.plot('radiance_noslit', wunits='cm-1', Iunits='W/m2/sr/cm_1')
        s.apply_slit(5)
        s.plot('radiance')
        w, t = s.get('transmittance_noslit')  # for use in multi-slabs configs

    Any tuple of numpy arrays (w, I) can also be converted into a Spectrum object
    from the :class:`~radis.spectrum.spectrum.Spectrum` class directly, or using 
    the :func:`~radis.spectrum.spectrum.Spectrum.calculated_spectrum` function. 
    All the following methods are equivalent::

        from radis import Spectrum, calculated_spectrum
        s1 = calculated_spectrum(w, I, wunit='nm', Iunit='mW/cm2/sr/nm')
        s2 = Spectrum.from_array(w, I, 'radiance_noslit', 
                               waveunit='nm', unit='mW/cm2/sr/nm')
        s3 = Spectrum({'radiance_noslit': (w, I)}, 
                      units={'radiance_noslit':'mW/cm2/sr/nm'},
                      waveunit='nm')

    Spectrum objects can be stored, retrieved, rescaled, resampled::

        from radis import load_spec
        s = load_spec('co_calculation.spec')
        s.rescale_path_length(0.5)                  # calculate for new path_length
        s.rescale_mole_fraction(0.02)   # calculate for new mole fraction
        s.resample(w_new)               # resample on new wavespace
        s.store('co_calculation2.spec')


    Notes
    -----
    
    Implementation:
    
        quantities are stored in self._q and self._q_conv dictionaries. They are better accessed
        with the get() function that deals with units and wavespace
        Note: we may move to a unique storage under self._q: in that case q['wavespace']
        and q_conv['wavespace'] need different names: ex. 'wavespace' and 'wavespace_conv'

    Wavebase:
    
        quantites are stored either in wavenum or wavelength base, but this doesnt
        matter as they are retrieved / plotted with get() and plot() functions
        which have units as input arguments


    Attributes
    ----------
    
    conditions : dict
        Stores computation / measurement conditions
        
    populations: dict
        Stores molecules, isotopes, electronic states and vibrational or 
        rovibrational populations
    
    
    See Also
    --------
    
    :func:`~radis.spectrum.spectrum.calculated_spectrum`, 
    :func:`~radis.spectrum.spectrum.transmittance_spectrum`, 
    :func:`~radis.spectrum.spectrum.experimental_spectrum`
    :meth:`~radis.spectrum.spectrum.from_array`
    :meth:`~radis.spectrum.spectrum.from_txt`
    :func:`~radis.tools.database.load_spec`
    
    
    '''

    # hardcode attribute names, but can save a lot of memory if hundreds of spectra
    __slots__ = ['_q', '_q_conv', 'units', 'conditions', 'cond_units', 'populations',
                  'lines', 'name', '_slit']

    def __init__(self, quantities, units=None, conditions=None, cond_units=None,
                 populations=None, lines=None, waveunit=None, wavespace=None, # deprecated
                 name=None, warnings=True):

        # Check inputs
        if wavespace is not None:
            warn(DeprecationWarning('wavespace replaced with waveunit'))
            waveunit = wavespace
        # ... Replace None attributes with dictionaries
        if conditions is None:
            conditions = {}
        if units is None:
            units = {}
        if cond_units is None:
            cond_units = {}
        if populations is None:
            populations = {}
        self._init_annotations()       # typing hints for get() and plot()

        # Deal with deprecated inputs
        # ... wavespace renamed waveunit
        if 'wavespace' in conditions:
            warn(DeprecationWarning('wavespace key in conditions is now named waveunit'))
            conditions['waveunit'] = conditions['wavespace']
            del conditions['wavespace']

        # Waveunit
        # ... Cast in standard waveunit format
        if waveunit is None and not 'waveunit' in conditions:
            raise AssertionError("waveunit ('nm', 'cm-1'?) has to be defined in `conditions`"+\
                                 "or with waveunit=")
        if waveunit is not None:
            waveunit = cast_waveunit(waveunit)
        if 'waveunit' in conditions:
            conditions['waveunit'] = cast_waveunit(conditions['waveunit'])
        # ... Make sure unit match
        if 'waveunit' in conditions:
            if waveunit is not None and conditions['waveunit'] != waveunit:
                raise ValueError('waveunit defined in conditions ({0}) and directly ({1}) dont match'.format(
                                            conditions['waveunit'], waveunit))
        else: # ... or define them in dictionary
            conditions['waveunit'] = waveunit

        # Check quantities format
        if len(quantities) == 0:
            raise AssertionError('Spectrum is created with no quantities. Add '+\
                                 '`radiance`, `transmittance` etc... with dict '+\
                                 'format. e.g: {`radiance`: (w,I)}')

        self._q = {}                # non convoluted quantities
        self._q_conv = {}           # convoluted quantities
        
        self._slit = {}            # hold slit function

        for k, v in quantities.items():
            try:
                assert len(v)==2
            except AssertionError:
                raise AssertionError('Attributes should have format (wavenumber, '\
                                     +'quantity) . Error with `{0}`'.format(k))

            w, I = v

            self._add_quantity(k, w, I, warnings=warnings)        # creates a copy of w,I

        # Finally, add our attributes
        self.conditions = conditions
        self.populations = populations
        self.lines = lines
        self.units = units
        self.cond_units = cond_units
        self.name = name

    # %% Constructors
    
    @classmethod
    def from_array(cls, w, I, quantity, waveunit, unit, *args, **kwargs):
        """
        Construct Spectrum from 2 arrays

        Parameters
        ----------
        
        w, I: array
            waverange and vector 
        
        quantity: str
            spectral quantity name
            
        waveunit: 'nm', 'cm-1'
            unit of waverange
            
        unit: str
            spectral quantity unit (arbitrary). Ex: 'mW/cm2/sr/nm' for radiance_noslit
            
        *args, **kwargs
            see :class:`~radis.spectrum.spectrum.Spectrum` doc
        
        Other Parameters
        ----------------
                
        Optional parameters:
    
        conditions: dict
            physical conditions and calculation parameters
    
        cond_units: dict
            units for conditions
    
        populations: dict
            a dictionary of all species, and levels. Should be compatible with other
            radiative codes such as Specair output. Suggested format:
            {molecules: {isotopes: {elec state: rovib levels}}}
            e.g::
            
                {'CO2':{1: 'X': df}}   # with df a Pandas Dataframe
            
    
        lines: pandas Dataframe
            all lines in databank (necessary for LineSurvey). Warning if you want to
            play with the lins content: The signification of columns in `lines` may be
            specific to a database format. Plus, some additional columns may have been
            added by the calculation (e.g: `Ei` and `S` for emission integral and
            linestrength in SpectrumFactory). Refer to the code to know what they mean
            (and their units)

        Returns
        -------
        :class:`~radis.spectrum.spectrum.Spectrum` object
        
        
        Examples
        --------
        
        Create a spectrum::
        
            from radis import Spectrum
            s = Spectrum.from_array(w, I, 'radiance_noslit', 
                                   waveunit='nm', unit='mW/cm2/sr/nm')
        
        
        See Also
        --------
        
        :func:`~radis.spectrum.spectrum.calculated_spectrum`, 
        :func:`~radis.spectrum.spectrum.transmittance_spectrum`, 
        :func:`~radis.spectrum.spectrum.experimental_spectrum`
        :meth:`~radis.spectrum.spectrum.from_txt`
        :func:`~radis.tools.database.load_spec`
        """
        
        # TODO: if medium not defined and quantities are given in cm-1, use vacuum
        
        quantities = {quantity:(w,I)}
        units = {quantity:unit}
        
        return cls(quantities, units, waveunit=waveunit, *args, **kwargs)

    @classmethod
    def from_txt(cls, file, quantity, waveunit, unit, *args, **kwargs):
        """
        Construct Spectrum from txt file
        

        Parameters
        ----------
        
        file: str
        
        quantity: str
            spectral quantity name
            
        unit: str
            spectral quantity unit
            
        *args, **kwargs
            the following inputs are forwarded to loadtxt: 'delimiter', 'skiprows'
            The rest if forwarded to Spectrum. see :class:`~radis.spectrum.spectrum.Spectrum`
            doc
        
        
        Other Parameters
        ----------------
                
        Optional parameters:
    
        conditions: dict
            physical conditions and calculation parameters
    
        cond_units: dict
            units for conditions
    
        populations: dict
            a dictionary of all species, and levels. Should be compatible with other
            radiative codes such as Specair output. Suggested format:
            {molecules: {isotopes: {elec state: rovib levels}}}
            e.g::
            
                {'CO2':{1: 'X': df}}   # with df a Pandas Dataframe
            
    
        lines: pandas Dataframe
            all lines in databank (necessary for LineSurvey). Warning if you want to
            play with the lins content: The signification of columns in `lines` may be
            specific to a database format. Plus, some additional columns may have been
            added by the calculation (e.g: `Ei` and `S` for emission integral and
            linestrength in SpectrumFactory). Refer to the code to know what they mean
            (and their units)


        Returns
        -------
        :class:`~radis.spectrum.spectrum.Spectrum` object
        
        
        Examples
        --------
        
        Generate an experimental spectrum from txt::
        
            from radis import Spectrum
            s = Spectrum.from_txt('spectrum.csv', 'radiance', waveunit='nm',
                                      unit='W/cm2/sr/nm')
        
        
        Notes
        -----
        
        Internally, the numpy ``loadtxt`` function is used and transposed::
        
            w, I = np.loadtxt(file).T

        
        See Also
        --------
        
        :func:`~radis.spectrum.spectrum.calculated_spectrum`, 
        :func:`~radis.spectrum.spectrum.transmittance_spectrum`, 
        :func:`~radis.spectrum.spectrum.experimental_spectrum`
        :meth:`~radis.spectrum.spectrum.from_array`
        :func:`~radis.tools.database.load_spec`
        
        """
        
        # Get input for loadtxt
        kwloadtxt = {}
        for k in ['delimiter', 'skiprows']:
            if k in kwargs:
                kwloadtxt[k] = kwargs.pop(k)
        
        w, I = np.loadtxt(file, **kwloadtxt).T
        quantities = {quantity:(w,I)}
        units = {quantity:unit}
        
        return cls(quantities, units, waveunit=waveunit, *args, **kwargs)


    # Public functions
    # %% ======================================================================
    # ----------------
    # XXX =====================================================================


    def get(self, var, wunit='nm', Iunit='default', medium='default',
            xunit=None, yunit=None,  # deprecated for wunit, Iunit
            copy=True):
        ''' Retrieve a spectral quantity from a Spectrum object. You can change 
        wavespace unit, intensity unit, or propagation medium.
        
        
        Parameters    
        ----------

        var: variable ('absorbance', 'transmittance', etc.) For full list see
            get_vars()

        wunit: 'cm', 'nm'
            wavelength / wavenumber unit. Default nm

        Iunit: unit for variable `var`
            if 'default', default unit for quantity `var` is used. See Spectrum.units
            to get the units. for radiance, one can use per wavelength (~ 'W/m2/sr/nm')
            or per wavenumber (~ 'W/m2/sr/cm_1') units

        medium: 'air', 'vacuum', 'default'
            returns wavelength as seen in air, or vacuum. If 'default' the
            value set in conditions is used. If None you better be sure of
            what you're doing.

        Other Parameters
        ----------------
        
        copy: boolean
            if True, returns a copy of the stored quantity (modifying it wont
            change the Spectrum object). Default True.
        
        Todo
        ----
        TODO: allow to get radiance in (W/sr/cm2/nm) or (W/sr/cm2) multiplying
        by FWHM.

    
        Returns
        -------

        w, I: array-like
            wavespace, quantity (ex: wavelength, radiance_noslit). For numpy
            users, note that these are copies (values) of the Spectrum quantity
            and not a view (reference): if you modify them the Spectrum is not
            changed

        '''

        # deprecation
        if xunit is not None:
            warn(DeprecationWarning('xunit replaced by wunit in Spectrum class'))
            wunit = xunit
        if yunit is not None:
            warn(DeprecationWarning('yunit replaced by Iunit in Spectrum class'))
            Iunit = yunit

        # check input
        if not var in self.get_vars():
            if var+'_noslit' in self.get_vars():
                raise KeyError(('`{0}` doesnt exist but `{1}` does'.format(var,
                                      var+'_noslit')) + '. Have you used .apply_slit()?')
            else:
                raise KeyError('{0} not in quantity list: {1}'.format(var, self.get_vars()))

        # Get quantity
        if var in CONVOLUTED_QUANTITIES:
            vartype = 'convoluted'
            I = self._q_conv[var]
        else:  # non convoluted
            vartype = 'non_convoluted'
            I = self._q[var]
            
        # Copy
        if copy: I = I.copy()
            
        # Get wavespace (in correct unit, and correct medium)
        wunit = cast_waveunit(wunit)
        if wunit == 'cm-1':
            w = self.get_wavenumber(vartype, copy=copy)
        else: # wunit == 'nm':
            w = self.get_wavelength(medium, vartype, copy=copy)
            
        # Convert y unit if necessary
        Iunit0 = self.units[var]
        if Iunit != 'default' and Iunit != Iunit0:
            if var in ['radiance', 'radiance_noslit']:
                # deal with the case where we want to get a radiance in per
                # wavelength unit (~ W/sr/cm2/nm) in wavenumber units (~ W/sr/cm2/cm_1),
                # or the other way round
                I = convert_universal(I, Iunit0, Iunit, self,
                                              per_nm_is_like='mW/sr/cm2/nm',
                                              per_cm_is_like='mW/sr/cm2/cm_1')
            elif var in ['emisscoeff']:
                # idem for emisscoeff in (~ W/sr/cm3/nm) or (~ /sr/cm3/cm_1)
                I = convert_universal(I, Iunit0, Iunit, self,
                                          per_nm_is_like='mW/sr/cm3/nm',
                                          per_cm_is_like='mW/sr/cm3/cm_1')
            elif var in ['absorbance']:  # no unit
                assert Iunit in ['', '-ln(I/I0)']
                # dont change the variable: I has no dimension
            elif var in ['transmittance']:  # no unit
                assert Iunit in ['', 'I/I0']
                # dont change the variable: I has no dimension
            else:
                I = conv2(I, Iunit0, Iunit)
        else:
            pass

        return w, I

    def _get_wavespace(self, which='any', copy=True):
        ''' Return wavespace (if the same for all quantities)

        
        Parameters    
        ----------

        which: 'convoluted', 'non_convoluted', 'any'
            return wavelength for convoluted quantities, non convoluted quantities, 
            or any. If any and both are defined, they have to be the same else 
            an error is raised. Default any.

        Other Parameters
        ----------------
        
        copy: boolean
            if True, returns a copy of the stored waverange (modifying it wont
            change the Spectrum object). Default True.
        
    
        Returns
        -------

        w: array_like
            (a copy of) spectrum wavespace for convoluted or non convoluted
            quantities

        '''

        if which == 'any':
            q_defined = 'wavespace' in self._q
            q_conv_defined = 'wavespace' in self._q_conv
            if q_defined and q_conv_defined:
                if not len(self._q['wavespace']) == len(self._q_conv['wavespace']):
                    raise ValueError('All wavespace not equal for calculated '+\
                                     "quantities. Can't use get_wavespace()")
                if not np.allclose(self._q['wavespace'], self._q_conv['wavespace']):
                    raise ValueError('All wavespace not equal for calculated '+\
                                     "quantities. Can't use get_wavespace()")
                w = self._q['wavespace']
            elif q_defined:
                w = self._q['wavespace']
            else:
                w = self._q_conv['wavespace']

        elif which == 'convoluted':
            w = self._q_conv['wavespace']

        elif which == 'non_convoluted':
            w = self._q['wavespace']

        else:
            raise ValueError("which has to be one of 'convoluted', 'non_convoluted', "+\
                             "'any'. Got {0}".format(which))

        if copy: w = w.copy()
        
        return w

    def get_wavelength(self, medium='default', which='any', copy=True):
        ''' Return wavelength in defined medium 
        
    
        Parameters    
        ----------

        which: 'convoluted', 'non_convoluted', 'any'
            return wavelength for convoluted quantities, non convoluted quantities, 
            or any. If any and both are defined, they have to be the same else 
            an error is raised. Default any.

        medium: 'air', 'vacuum', 'default'
            returns wavelength as seen in air, or vacuum. If 'default' the
            value set in conditions is used. If None you better be sure of
            what you're doing.

        Other Parameters
        ----------------
        
        copy: boolean
            if True, returns a copy of the stored waverange (modifying it wont
            change the Spectrum object). Default True.
        
    

        Returns
        -------

        w: array_like
            (a copy of) spectrum wavelength for convoluted or non convoluted
            quantities
        '''
        
        # Check input        
        if not medium in ['default', 'air', 'vacuum']:
            raise NotImplementedError('Unknown propagating medium: {0}'.format(medium))
        stored_medium = self.conditions.get('medium', None)
        if not stored_medium in ['air', 'vacuum', None]:
            raise NotImplementedError('Unknown propagating medium: {0}'.format(stored_medium))
            
        # Now get wavespace
        w = self._get_wavespace(which=which, copy=copy)
        if self.get_waveunit() == 'cm-1':
            w = cm2nm(w)       # we get vacuum wavelength
        
            # Correct for propagation medium (air, vacuum)
            if medium == 'default':
                if stored_medium is None:
                    raise ValueError('Medium not defined in Spectrum conditions. We cant '+\
                                     'guess in which medium to show wavelengths. '+\
                                     "Update conditions, explicitely use medium="+\
                                     "'air'/'vacuum', or work with wavenumbers "+\
                                     "(wunit='cm-1')")
                else:
                    medium = stored_medium
            
            # Now medium is either 'air' or 'vacuum'
            if medium == 'air':
                w = vacuum2air(w)  
            else:   # medium == 'vacuum'
                pass  # no change needed

        else: # self.get_waveunit() == 'nm'
            # we got w in 'stored_medium'
            
            # Correct for propagation medium (air, vacuum) depending on the value of stored_medium
            if stored_medium == 'vacuum':
                if medium == 'air':
                    w = vacuum2air(w)       # w is in nm since `w = cm2nm(w_cm)`
                else:
                    pass # no change needed
            elif stored_medium == 'air':
                if medium == 'vacuum':
                    w = air2vacuum(w)
                else:
                    pass   # no change needed
            else:  # stored_medium is None:
                if medium != 'default':
                    raise ValueError('Medium not defined in Spectrum conditions. We cant '+\
                                     'convert wavelengths to {0}. '.format(medium)+\
                                     "Update conditions, or use medium='default'")
                else:
                    pass # no change needed
            
        return w

    def get_wavenumber(self, which='any', copy=True):
        ''' Return wavenumber (if the same for all quantities)
        
    
        Parameters    
        ----------

        which: 'convoluted', 'non_convoluted', 'any'
            return wavenumber for convoluted quantities, non convoluted quantities, 
            or any. If any and both are defined, they have to be the same else 
            an error is raised. Default any.


        Other Parameters
        ----------------
        
        copy: boolean
            if True, returns a copy of the stored waverange (modifying it wont
            change the Spectrum object). Default True.
        
    
        Returns
        -------

        w: array_like
            (a copy of) spectrum wavenumber for convoluted or non convoluted
            quantities
        '''
        w = self._get_wavespace(which=which, copy=copy)
        
        if self.get_waveunit() == 'nm':     # we want w in vacuum before converting
            # Correct for propagation medium (air, vacuum)
            stored_medium = self.conditions.get('medium', None)
            if stored_medium is None:
                raise KeyError('Medium not defined in Spectrum conditions. We cant '+\
                                 'derive wavenumber if we dont know whether '+\
                                 'wavelengths are in vacuum or air. '+\
                                 "Update conditions")
            elif stored_medium == 'vacuum':
                pass
            elif stored_medium == 'air':
                w = air2vacuum(w)
            else:
                raise NotImplementedError('Unknown propagating medium: {0}'.format(stored_medium))

            # Convert to wavenumber
            w = nm2cm(w)
        return w

    def get_radiance(self, Iunit='mW/cm2/sr/nm', copy=True):
        ''' Return radiance in whatever unit, and can even convert from ~1/nm
        to ~1/cm-1 (and the other way round)
        
        
        Other Parameters
        ----------------
        
        copy: boolean
            if True, returns a copy of the stored waverange (modifying it wont
            change the Spectrum object). Default True.
        
    
        '''
        
        return self.get('radiance', Iunit=Iunit, copy=copy)[1]

#        Iunit0 = self.units['radiance']
#        
#        w_cm, I = self.get('radiance', wunit='cm-1', Iunit=Iunit0, copy=copy)
#
#        return convert_universal(I, Iunit0, Iunit, w_cm,
#                                     per_nm_is_like='mW/sr/cm2/nm',
#                                     per_cm_is_like='mW/sr/cm2/cm_1')

    def get_radiance_noslit(self, Iunit='mW/cm2/sr/nm', copy=True):
        ''' Return radiance (non convoluted) in whatever unit, and can even
        convert from ~1/nm to ~1/cm-1 (and the other way round)
        
        
        Other Parameters
        ----------------
        
        copy: boolean
            if True, returns a copy of the stored waverange (modifying it wont
            change the Spectrum object). Default True.
        
    
        '''
        
        return self.get('radiance_noslit', Iunit=Iunit, copy=copy)[1]

#        Iunit0 = self.units['radiance_noslit']
#        w_cm, I = self.get('radiance_noslit', wunit='cm-1', Iunit=Iunit0, copy=copy)
#
#        return convert_universal(I, Iunit0, Iunit, w_cm,
#                                     per_nm_is_like='mW/sr/cm2/nm',
#                                     per_cm_is_like='mW/sr/cm2/cm_1')

    def get_name(self):
        ''' Return Spectrum name. If not defined, returns 'spectrum{id}' with
        the Python id object
        '''
        try:
            self.name
        except AttributeError:
            warn('Spectrum has no .name attribute and is probably outdated. Update!')
            self.name = None

        if self.name is not None:
            name = self.name
        else:
            name = 'spectrum{0}'.format(id(self))

        return name
    
    def savetxt(self, filename, var, wunit='nm', Iunit='default', medium='default'):
        ''' Export spectral quantity var to filename
        
        Parameters    
        ----------
        
        filename: str
            file name
        
        var: str
            which spectral variable ot export 
            
        Other Parameters
        ----------------
        
        wunit, Iunit, medium: str
            see :meth:`~radis.spectrum.spectrum.Spectrum.get` for more information
        
        Notes
        -----
        
        Export variable as::
            
            np.savetxt(filename, np.vstack(self.get(var, wunit=wunit, Iunit=Iunit,
                                                medium=medium)).T, header=header)
        
        See Also
        --------
        
        :meth:`~radis.spectrum.spectrum.Spectrum.store`
        
        '''
        
        # Get units to export 
        wunit = cast_waveunit(wunit)
        wmedium = ' [{0}]'.format(medium) if medium != 'default' else ''
        if wunit == 'cm-1':
            xlabel='Wavenumber (cm-1)'
        else: # wunit == 'nm'
            xlabel='Wavelength (nm){0}'.format(wmedium)

        if Iunit == 'default':
            try:
                yunit = self.units[var]
            except KeyError:  # unit not defined in dictionary
                yunit = 'a.u'
        else:
            yunit = Iunit
            
        header = '{0}\t{1} ({2})'.format(xlabel, var, yunit)

        np.savetxt(filename, np.vstack(self.get(var, wunit=wunit, Iunit=Iunit,
                                                medium=medium)).T, header=header)
    
    def update(self, quantity='all', optically_thin='default', verbose=True):
        ''' Calculate missing quantities: ex: if path_length and emisscoeff
        are given, recalculate radiance_noslit
        
    
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

        return update(self, quantity=quantity, optically_thin=optically_thin,
                      verbose=verbose)

    # Rescale functions

    def rescale(self, new_path_length, old_path_length=None,
                force=False):
        ''' Deprecated. See rescale_path_length '''

        warn(DeprecationWarning("rescale have been replaced by more explicit "+\
                                "rescale_path_length"))

        return self.rescale_path_length(new_path_length, old_path_length=old_path_length,
                force=force)

    def rescale_path_length(self, new_path_length, old_path_length=None,
                force=False):
        ''' Rescale spectrum to new path length. Starts from absorption coefficient
        and emission coefficient, and solves the RTE again for the new path length
        Convoluted values (with slit) are dropped in the process.

    
        Parameters    
        ----------

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
        
        return rescale_path_length(self, new_path_length=new_path_length, 
                                   old_path_length=old_path_length,
                                   force=force)

    def rescale_mole_fraction(self, new_mole_fraction, old_mole_fraction=None,
                ignore_warnings=False, force=False):
        ''' Update spectrum with new molar fraction
        Convoluted values (with slit) are dropped in the process.

    
        Parameters    
        ----------

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
        ----

        Add warning when too large rescaling
        '''
        
        return rescale_mole_fraction(self, new_mole_fraction=new_mole_fraction, 
                                     old_mole_fraction=old_mole_fraction,
                                     ignore_warnings=ignore_warnings, 
                                     force=force)

    def get_integral(self, var, wunit='nm', Iunit='default', **kwargs):
        ''' Returns integral of variable 'var' '''
        
        w, I = self.get(var, wunit=wunit, Iunit=Iunit, **kwargs)
        return abs(np.trapz(I, x=w))

    def power(self, unit='mW/cm2/sr'):
        ''' Returns integrated radiance power density'''
        
        raise DeprecationWarning('Spectrum.power() deprecated. Use get_power() instead')

    def get_power(self, unit='mW/cm2/sr'):
        ''' Returns integrated radiance (no slit) power density'''

        P = self.get_integral('radiance_noslit', wunit='nm', Iunit='mW/cm2/sr/nm')
        # P is in mW/cm2/sr/nm * nm
        return conv2(P, 'mW/cm2/sr', unit)

    # %% Plotting routines

    def get_vars(self, which='any'):
        ''' Returns all spectral quantities stored in this object (convoluted or 
        non convoluted)
        
        Parameters
        ----------
        
        which: 'any', 'convoluted', 'non convoluted' 
        
        '''
        if which == 'any':
            varlist = list(self._q.keys())+list(self._q_conv.keys())
        elif which == 'convoluted':
            varlist = list(self._q_conv.keys())
        elif which == 'non convoluted':
            varlist = list(self._q.keys())
        else:
            raise ValueError('Unexpected value: {0}'.format(which))
            
        # remove wavespace
        varlist = [k for k in varlist if k != 'wavespace']
        return varlist

    def _get_items(self):
        ''' Return a dictionary of tuples, e.g, {'radiance':(w,I), 'transmittance_noslit':(w_ns,T)}
        
        In the general case, users should use s.get()
        '''
        items = {k:(self._q['wavespace'],v) for k,v in self._q.items() if k != 'wavespace'}
        items.update({k:(self._q_conv['wavespace'],v) for k,v in self._q_conv.items() if k != 'wavespace'})
        return items

    def plot(self, var=None, wunit='default', Iunit='default', show_points=False,
             nfig=None, yscale='linear', medium='default',
             xunit=None, yunit=None, # deprecated
             normalize=False,
             **kwargs):
        '''
    
        Parameters    
        ----------

        var: variable (`absorbance`, `transmittance`, `transmittance_noslit`, etc.)
            For full list see get_vars(). If None, plot the first thing in
            the Spectrum (with a preference for non convoluted quantities)

        wunit: `default`, `cm-1`, `nm`
            wavelength or wavenumber unit. If `default`, Spectrum waveunit is used.

        Iunit: unit for variable
            if `default`, default unit for quantity `var` is used.
            for radiance, one can use per wavelength (~ `W/m2/sr/nm`) or
            per wavenumber (~ `W/m2/sr/cm_1`) units


        Other Parameters
        ----------------
        
        Plot parameters inputs:

        show_points: boolean
            show calculated points. Default True

        nfig: int, None, or 'same'
            plot on a particular figure. 'same' plots on current figure.

        yscale: 'linear', 'log'
            plot yscale

        normalize: boolean
            option to normalize quantity to 1 (ex: for radiance). Default False

        **kwargs: **dict
            kwargs forwarded as argument to plot (e.g: lineshape
            attributes: `lw=3, color='r'`)

        '''

        # Check inputs, get defaults
        # ------
        if xunit is not None:
            warn(DeprecationWarning('xunit replaced with wunit'))
            wunit = xunit
        if yunit is not None:
            warn(DeprecationWarning('yunit replaced with Iunit'))
            Iunit = yunit
        if var in ['intensity', 'intensity_noslit']:
            raise ValueError('`intensity` not defined. Use `radiance` instead')

        if var is None:    # if nothing is defined, try these first:
            params = self.get_vars()
            if 'radiance' in params:
                var = 'radiance'
            elif 'radiance_noslit' in params:
                var = 'radiance_noslit'
            elif 'transmittance' in params:
                var = 'transmittance'
            elif 'transmittance_noslit' in params:
                var = 'transmittance_noslit'
            else:
                # or plot the first variable we find
                var = list(params)[0]
                if var.replace('_noslit', '') in params:
                    var = var.replace('_noslit', '')

        if wunit == 'default':
            wunit = self.get_waveunit()
        wunit = cast_waveunit(wunit)
            
        # Get variable
        x, y = self.get(var, wunit=wunit, Iunit=Iunit, medium=medium)
        
        wmedium = ' [{0}]'.format(medium) if medium != 'default' else ''
        if wunit == 'cm-1':
            xlabel='Wavenumber (cm-1)'
        else: # wunit == 'nm'
            xlabel='Wavelength (nm){0}'.format(wmedium)

        if Iunit == 'default':
            try:
                Iunit0 = self.units[var]
            except KeyError:  # unit not defined in dictionary
                Iunit0 = 'a.u'
            Iunit = Iunit0
            
        # cosmetic changes
        Iunit = make_up(Iunit)

        # Plot
        # -------
        if normalize:
            y /= y.max()
            Iunit = 'norm'
            
        set_style('origin')
        if nfig == 'same':
            nfig = plt.gcf().number
        plt.figure(nfig)
        # Add extra plotting parameters
        if 'lw' not in kwargs and 'linewidth' not in kwargs:
            kwargs['lw'] = 0.5
        line, = plt.plot(x, y, **kwargs)  # note: '-k' by default with style origin for first plot
        if show_points: plt.plot(x, y,'o', color='lightgrey', **kwargs)
        plt.ticklabel_format(useOffset=False, axis='x')
        plt.xlabel(xlabel)
        plt.ylabel(make_up('{0} ({1})'.format(var, Iunit)))

        plt.yscale(yscale)

        if 'label' in kwargs:
            plt.legend()

        fix_style(str('origin'))

        plt.show()

        return line
    
    def get_populations(self, molecule=None, isotope=None, electronic_state=None):
        ''' Return populations that are featured in the spectrum, either as 
        upper or lower levels 
        
    
        Parameters    
        ----------
        
        molecule: str, or None
            if None, only one molecule must be defined. Else, an error is raised
            
        isotope: int, or None
            isotope number. if None, only one isotope must be defined. Else, 
            an error is raised
            
        electronic_state: str
            if None, only one electronic state must be defined. Else, an error 
            is raised
            

        Returns
        -------
        
        pandas dataframe of levels, where levels are the index, 
        and 'Evib' and 'nvib' are featured 
        
        Notes
        -----
        
        Structure:
            
            {molecule: {isotope: {electronic_state: {'vib': pandas Dataframe,    # (copy of) vib levels
                                                     'rovib': pandas Dataframe,  # (copy of) rovib levels
                                                     'Ia': float    # isotopic abundance
                                                     }}}}
            
        (If Spectrum generated with RADIS, structure should match that of 
        SpectrumFactory.get_populations())
        
            
        '''
        
        # Check inputs, get default values
        populations = self.populations
        if populations is None or populations == {}:
            raise ValueError('Populations not defined')
        if type(populations) != dict:
            raise TypeError('Method defined for populations as dictionary')
        if molecule is None:
            if len(list(populations.keys()))!=1:
                raise ValueError('Please choose which molecule among: {0}'.format(
                        list(populations.keys())))
            molecule = list(populations.keys())[0]
        if isotope is None:
            if len(list(populations[molecule].keys()))!=1:
                raise ValueError('Please choose which isotope among: {0}'.format(
                        list(populations[molecule].keys())))
            isotope = list(populations[molecule].keys())[0]
        if electronic_state is None:
            if len(list(populations[molecule][isotope].keys()))!=1:
                raise ValueError('Please choose which electronic state among: {0}'.format(
                        list(populations[molecule][isotope].keys())))
            electronic_state = list(populations[molecule][isotope].keys())[0]
            
        if __debug__: printdbg('get vib populations for {0}({1})[iso{2}]'.format(
                molecule, electronic_state, isotope))
        
        # Return            
        return populations[molecule][isotope][electronic_state]
        
    def get_vib_levels(self, molecule=None, isotope=None, electronic_state=None,
                       first=None):
        ''' Return vibrational levels in the spectrum (energies, populations)
        
    
        Parameters  
        ----------  
        
        molecule: str, or None
            if None, only one molecule must be defined. Else, an error is raised
            
        isotope: int, or None
            isotope number. if None, only one isotope must be defined. Else, 
            an error is raised
            
        electronic_state: str
            if None, only one electronic state must be defined. Else, an error 
            is raised
            
        first: int, or 'all' or None
            only show the first N levels. If None or 'all', all levels are shown
            

        Returns
        -------
        
        pandas dataframe of levels, where levels are the index, 
        and 'Evib' and 'nvib' are featured 
            
        '''
        
        pops = self.get_populations(molecule=molecule, isotope=isotope,
                                     electronic_state=electronic_state)
        
        try:
            vib_pops = pops['vib']
        except KeyError:
            raise KeyError('Vibrational levels not defined in this Spectrum object. '+\
                           "If using RADIS, make sure you chose export_populations='vib'")
        
        if first is not None:
            if not 'nvib' in vib_pops:
                raise KeyError('Vibrational populations (nvib) not calculated in this '+\
                               'Spectrum object. Cant get first most populated levels. '+\
                               'If using RADIS, make sure you used a non_eq_spectrum '+\
                               'calculation')
            if first == 'all': first = None
            out = vib_pops.sort_values(by='nvib', ascending=False)[:first]
            
        else:
            out = vib_pops
        
        # Return            
        return out
        
    def get_rovib_levels(self, molecule=None, isotope=None, electronic_state=None,
                       first=None):
        ''' Return rovibrational levels calculated in the spectrum (energies, 
        populations)
        
    
        Parameters    
        ----------
        
        molecule: str, or None
            if None, only one molecule must be defined. Else, an error is raised
            
        isotope: int, or None
            isotope number. if None, only one isotope must be defined. Else, 
            an error is raised
            
        electronic_state: str
            if None, only one electronic state must be defined. Else, an error 
            is raised
            
        first: int, or 'all' or None
            only show the first N levels. If None or 'all', all levels are shown


        Returns
        -------
    
        pandas dataframe of levels, where levels are the index, 
        and 'Evib' and 'nvib' are featured 
            
        '''
        
        pops = self.get_populations(molecule=molecule, isotope=isotope,
                                     electronic_state=electronic_state)
        
        try:
            rovib_pops = pops['rovib']
        except KeyError:
            raise KeyError('Vibrational levels not defined in this Spectrum object. '+\
                           "If using RADIS, make sure you chose export_populations='vib'")
        
        if first is not None:
            if not 'n' in rovib_pops:
                raise KeyError('Rovibrational populations (n) not calculated in this '+\
                               'Spectrum object. Cant get first most populated levels. '+\
                               'If using RADIS, make sure you used a non_eq_spectrum '+\
                               'calculation')
            if first == 'all': first = None
            out = rovib_pops.sort_values(by='n', ascending=False)[:first]
            
        else:
            out = rovib_pops
        
        # Return            
        return out
        
    def plot_populations(self, what=None, nunit='', correct_for_abundance=False, **kwargs):
        ''' Plots vib populations if given and format is valid
        
    
        Parameters    
        ----------
        
        what: 'vib', 'rovib', None
            if None plot everything
        
        nunit: '', 'cm-3'
            plot either in a fraction of vibrational levels, or a molecule 
            number in in cm-3
        
        correct_for_abundance: boolean
            if True, multiplies each population by the isotopic abundance 
            (as it is done during the calculation of emission integral)
        
        kwargs: **dict
            are forwarded to the plot
        '''
        
        # Check input, get defaults
        pops = self.populations
        if not isinstance(pops, dict):
            raise ValueError('Populations not defined in given Spectrum')
        
        # Plot function
        def _plot_elec_state(what, df, state_name, Ia, fig):
            if what == 'vib':
                E, n, g = df['Evib'], df['nvib'], df['gvib']
                ylabel = 'n/g$_{vib}$'
                title = 'Vibrational populations'
            elif what == 'rovib':
                E, n, g = df['E'], df['n'], df['gj']
                ylabel = 'n/(2J+1)'
                title = 'Rovibrational populations'
                
            if correct_for_abundance:
                n = n * Ia
    
            if nunit == 'cm-3':
                from radis.phys.constants import k_b
                try:    
                    P_mbar = self.conditions['pressure_mbar']  # mbar
                    T = self.conditions['Tgas']
                    mfrac = self.conditions['mole_fraction']
                except KeyError:
                    raise KeyError('P_mbar (pressure), T (Tgas) and n (mole_fraction) '+\
                                   'are needed to calculate total number density in (cm-3)')
                N = P_mbar*1e2/k_b/T*mfrac*1e-6
                n = n * N
                unitlabel = ' [cm-3]'
            elif nunit == '':
                unitlabel = ' [fraction]'
            else:
                raise ValueError('Unknown unit: {0}'.format(nunit))
            
            # Plot
            if fig is None:
                fig = plt.figure()
                plt.xlabel('Energy (cm-1)')
                plt.ylabel(ylabel+unitlabel)
                plt.yscale('log')
                plt.title(title)
            ax = fig.gca()
            ax.plot(E, n/g, 'o', label=state_name, **kwargs)
            
            return fig

        # Initialize figures, styles
        fig_vib = None
        fig_rovib = None
        set_style('origin')
        
        # Loop over all molecules, all isotopes, all electronic states
        # Note that the below works for both dict and pandas dataframe
        
        for molecule, isotopes in pops.items():
            for isotope, elec_states in isotopes.items():
                for elec_state, content in elec_states.items():
                    state_name = '{0}({1})(iso{2})'.format(molecule, elec_state, 
                                  isotope)

                    Ia = None                    
                    if correct_for_abundance:
                        if 'Ia' in list(content.keys()):
                            Ia = content['Ia']
                        else:
                            raise KeyError('Ia: isotopic abundance not defined in '+\
                                           'Spectrum populations')
                    
                    for k in content.keys():
                        
                        if k == 'rovib' and (what == k or what is None):
                            df = self.get_rovib_levels(molecule, isotope, elec_state,
                                                       first='all')  # sort by n + check is defined)
                            fig_rovib = _plot_elec_state('rovib', df, state_name, 
                                                         Ia, fig_rovib)
                            
                        if k == 'vib' and (what == k or what is None):
                            df = self.get_vib_levels(molecule, isotope, elec_state,
                                                     first='all')  # sort by n + check is defined)
                            fig_vib = _plot_elec_state('vib', df, state_name, 
                                                       Ia, fig_vib)
                            
        # Update
        for fig in [fig_vib, fig_rovib]:
            if fig is not None:
                ax = fig.gca()
                ax.legend()
                fix_style('origin', ax)

    # %% ------------------ Instrumental Slit Function ---------------------

    def apply_slit(self, slit_function, unit='nm', shape='triangular',
                   center_wavespace=None, plot_slit=False, norm_by='area',
                   mode='valid', store=True,
                   *args, **kwargs):
        ''' Apply a slit function to all quantities in Spectrum. Slit function
        can be generated with usual shapes (see `shape=`) or imported from an
        experimental slit function.

        Warning with units: read about `unit` and `return_unit` parameters.

    
        Parameters    
        ----------

        slit_function: float (nm) / str (.txt file)
            If float:
                generate slit function with FWHM of slit function (in nm or
                cm-1 depending on slit_unit)
            If .txt file:
                import experimental slit function: format must be 2-columns with
                wavelengths and intensity (doesn't have to be normalized)

        unit: 'nm' or 'cm-1'
            unit of slit_function (FWHM, or imported file)

        shape: 'triangular', 'trapezoidal', 'gaussian'
            which shape to use when generating a slit. Will call,
             respectively, :func:`~radis.tools.slit.triangular_slit`, 
             :func:`~radis.tools.slit.trapezoidal_slit`, 
             :func:`~radis.tools.slit.gaussian_slit`. Default 'triangular'

        center_wavespace: float, or None
            center of slit when generated (in unit). Not used if slit is imported.

        norm_by: 'area', 'max', 'max2'
            normalisation type:
            
            - 'area' normalizes the slit function to an area
              of 1. It conserves energy, and keeps the same units. 
            
            - 'max' normalizes the slit function to a maximum of 1. 
              The convoluted spectrum units change (they are 
              multiplied by the spectrum waveunit, e.g: a radiance 
              non convoluted in mW/cm2/sr/nm on a wavelength (nm).
              range will yield a convoluted radiance in mW/cm2/sr. 
              Note that the slit is set to 1 in the Spectrum wavespace
              (i.e: a Spectrum calculated in cm-1 will have a slit
              set to 1 in cm-1). 
            
            - 'max2': read about it in :func:`~radis.tools.slit.get_slit_function` 
              docstrings. 
            
            Default 'area'
    
        mode: 'valid', 'same'
            'same' returns output of same length as initial spectra, 
            but boundary effects are still visible. 'valid' returns 
            output of length len(spectra) - len(slit) + 1, for 
            which lines outside of the calculated range have
            no impact. Default 'valid'. 

        store: boolean
            if True, store slit in the Spectrum object so it can be retrieved with 
            :meth:`~radis.spectrum.spectrum.Spectrum.get_slit` and plot with 
            :meth:`~radis.spectrum.spectrum.Spectrum.plot_slit`. Default True

        Other Parameters
        ----------------

        *args, **kwargs
            are forwarded to slit generation or import function

        In particular:

        energy_threshold: float
             tolerance fraction when resampling. Default 1e-3 (0.1%)
             If areas before and after resampling differ by 
             more than that an error is raised. 


        Notes
        -----
        
        Implementation:

        Apply :func:`~radis.tools.slit.convolve_with_slit` 
        for all quantities in .vars that ends with
        _noslit. Generate a triangular instrumental slit function 
        (or any other shape depending of shape=) with base
        `slit_function_base` (Uses the central wavelength of the spectrum
        for the slit function generation)

        We deal with several special cases (which makes the code 
        a little heavy, but the method very versatile):
            
        - when slit unit and spectrum unit arent the same
        - when spectrum is not evenly spaced


        Examples
        --------

        To manually apply the slit to a particular quantity use::

            wavenum, quantity = s['quantity']
            s['convolved_quantity'] = convolve_slit(wavenum, quantity,
                slit_function_base)

        See convolve_slit for more details on Units and Normalization

        The slit is made considering the "center wavelength" which is
        the mean wavelength of the full spectrum you are applying it to.


        See Also
        --------
        
        :func:`~radis.tools.slit.get_slit_function`, 
        :func:`~radis.tools.slit.convolve_with_slit`

        '''
        # TODO: add warning if FWHM >= wstep(spectrum)/5    
        
        from radis.tools.slit import convolve_with_slit, get_slit_function, cast_waveunit

        # Check inputs
        # ---------
        if 'slit_unit' in kwargs:
            unit = kwargs.pop('slit_unit')
            warn(DeprecationWarning('slit_unit was renamed unit'))

        unit = cast_waveunit(unit)

        varlist = [k for k in self.get_vars() if k.endswith('_noslit')]

        if len(varlist) == 0:
            raise AssertionError('No variables to apply slit on. Variable names '+\
                                 'to be convolved should end with _noslit')

        # Forward relevant inputs to convolution instead of slit function generation
        kwargsconvolve = {}
        for kw in ['slit_dispersion', 'verbose']:
            if kw in kwargs:
                kwargsconvolve.update({kw:kwargs.pop(kw)})

        # For non evenyly distributed cases we take the minimum wstep among the
        # spectral range (and resample in convolve_with_slit)
        # Note: by construction all variables now have the same wavespace
        w = self._q['wavespace']  # non convoluted wavespace
        wstep = abs(diff(w)).min()
        assert wstep>0
        waveunit = self.get_waveunit()
        
        if __debug__: printdbg('apply_slit: {0} in {1}, center `{2}`{1}, applied in waveunit {3}'.format(
                slit_function, unit, center_wavespace, waveunit))
        
        if center_wavespace is None:
            # center_wavespace should be ~ unit
            center_wavespace = w[len(w)//2]  # w ~ waveunit
            if waveunit == 'cm-1' and unit == 'nm':
                center_wavespace = cm2nm(center_wavespace)     # wavenum > wavelen
            elif waveunit == 'nm' and unit == 'cm-1':
                center_wavespace = nm2cm(center_wavespace)     # wavelen > wavenum

        # Get slit once and for all (and convert the slit unit 
        # to the Spectrum `waveunit` if wavespaces are different)
        # -------
        wslit, Islit = get_slit_function(slit_function, unit=unit, norm_by = norm_by,
                                         shape=shape, center_wavespace = center_wavespace,
                                         return_unit = waveunit, wstep=wstep,
                                         plot=plot_slit, *args, **kwargs)

        # Apply to all variables
        # ---------
        for qns in varlist:
            q = qns[:-7]   # new name  (minus '_noslit')
            I = self._q[qns]

            # Convolve and store the output in a new variable name (quantity name minus `_noslit`)
            # Create if requireds
            w_conv, I_conv =  convolve_with_slit(w, I, wslit, Islit, norm_by = None, # already norm.
                                          mode=mode, **kwargsconvolve)
            self._q_conv['wavespace'] = w_conv
            self._q_conv[q] = I_conv

            # Get units
            if norm_by == 'area':
                self.units[q] = self.units[qns]
            elif norm_by == 'max':
                new_unit = '{0}*{1}'.format(self.units[qns], unit.replace('cm-1', 'cm_1'))
                                        # because it's like if we multiplied 
                                        # by slit FWHM in the wavespace it was
                                        # generated
                # simplify unit:
                try:
                    new_unit = '{:~P}'.format(Q_(new_unit).units)
                except UndefinedUnitError:
                    pass
                self.units[q] = new_unit
            elif norm_by == 'max2':   # difference with 'max': unit is multiplied by [unit] not [return_unit]
                new_unit = '{0}*{1}'.format(self.units[qns], waveunit.replace('cm-1', 'cm_1'))
                                        # because it's like if we multiplied 
                                        # by slit FWHM in the wavespace it was
                                        # generated
                # simplify unit:
                try:
                    new_unit = '{:~P}'.format(Q_(new_unit).units)
                except UndefinedUnitError:
                    pass
                self.units[q] = new_unit
            else:
                raise ValueError('Unknown normalization type: {0}'.format(norm_by))
                
        # Store slit in Spectrum, in the Spectrum unit
        if store:
            self._slit['wavespace'] = wslit  # in 'waveunit'
            self._slit['intensity'] = Islit

        # Update conditions
        self.conditions['slit_function'] = slit_function
        self.conditions['slit_unit'] = unit  # input slit unit
        self.conditions['norm_by'] = norm_by

        return
    
    def get_slit(self):
        ''' Get slit function that was applied to the Spectrum 
        
        Returns
        -------
        
        wslit, Islit: array
            slit function with wslit in waveunit. See 
            :meth:`~radis.spectrum.spectrum.Spectrum.get_waveunit`
            
        '''
    
        # Make sure that slit is stored already
        try:
            wslit = self._slit['wavespace'] # stored in Spectrum waveunit
            Islit = self._slit['intensity']
        except KeyError:
            raise KeyError('Slit function not found in Spectrum '+\
                             'conditions. Have you used Spectrum.apply_slit '+\
                             'with store=True?')
            
        return wslit, Islit
    
    def plot_slit(self, wunit=None):
        ''' Plot slit function that was applied to the Spectrum 
        
        Parameters
        ----------
        
        wunit: 'nm', 'cm-1', or None
            plot slit in wavelength or wavenumber. If None, use the unit
            the slit in which the slit function was given. Default None
        '''
        
        from radis.tools.slit import plot_slit
        
        # Check inputs
        assert wunit in ['nm', 'cm-1', None]
        if wunit is None:
            wunit = self.conditions['slit_unit']

        # Get slit arrays (in Spectrum.waveunit)
        wslit, Islit = self.get_slit()
        
        # Get slit unit
        norm_by = self.conditions['norm_by']
        waveunit = self.get_waveunit()
        if norm_by == 'area':
            Iunit = '1/{0}'.format(waveunit)
        elif norm_by == 'max': # set maximum to 1
            Iunit = '1'
        elif norm_by is None:
            Iunit=None
        else:
            raise ValueError('Unknown normalization type: `norm_by` = {0}'.format(norm_by))
        
        # Plot in correct unit  (plot_slit deals with the conversion if needed)
        return plot_slit(wslit, Islit, waveunit=waveunit, plot_unit=wunit, 
                         Iunit=Iunit)
        
        
    def line_survey(self, overlay=None, wunit='cm-1', cutoff=None, medium = 'default', 
                    xunit=None,  # deprecated
                    *args, **kwargs):
        ''' Plot Line Survey (all linestrengths used for calculation)
        Output in Plotly (html) 

    
        Parameters    
        ----------

        spec: Spectrum
            result from SpectrumFactory calculation (see spectrum.py)

        overlay: 'absorbance', 'transmittance', 'radiance', etc... or list of the above, or None
            overlay Linestrength with specified variable calculated in `spec`.
            Get the full list with the :meth:`~radis.spectrum.spectrum.Spectrum.get_vars`
            method. Default None. 

        wunit: 'nm', 'cm-1'
            wavelength / wavenumber units
            
        medium: {'air', 'vacuum', 'default'}
            Choose whether wavelength are shown in air or vacuum. If 'default'  
            lines are shown as stored in the spectrum. 

        Other inputs are passed to :func:`~radis.tools.line_survey.LineSurvey`

        Other Parameters
        ----------------
        
        ex:

        Iunit: `hitran`, `splot` 
            Linestrength output units:
                
            - `hitran`: (cm-1/(molecule/cm-2))
            - `splot`: (cm-1/atm)   (Spectraplot units [2]_)
            
            Note: if not None, cutoff criteria is applied in this unit.
            Not used if plot is not 'S'

        barwidth: float
            With of bars in LineSurvey. Default 0.07

        See :func:`~radis.tools.line_survey.LineSurvey` documentation

        
        Returns
        -------
        
        Plot in Plotly. See Output in [1]_
        
        
        Examples
        --------
            
        An example using the :class:`~radis.lbl.factory.SpectrumFactory` to generate a spectrum::
        
            from radis import SpectrumFactory
            sf = SpectrumFactory(
                                 wavenum_min=2380,
                                 wavenum_max=2400,
                                 mole_fraction=400e-6,
                                 path_length=100,  # cm
                                 isotope=[1],
                                 db_use_cached=True) 
            sf.load_databank('HITRAN-CO2-TEST')
            s = sf.eq_spectrum(Tgas=1500)
            s.apply_slit(0.5)
            s.line_survey(overlay='radiance_noslit', barwidth=0.01)
    
        
        References
        ----------
        
        .. [1] `RADIS Online Documentation (LineSurvey) <https://radis.readthedocs.io/en/latest/tools/line_survey.html>`__

        .. [2] `SpectraPlot <http://www.spectraplot.com/survey>`__


        See Also
        --------

        :func:`~radis.tools.line_survey.LineSurvey`
    
    
        '''
        
        from radis.tools.line_survey import LineSurvey

        # Check inputs
        if xunit is not None:
            warn(DeprecationWarning('xunit replaced with wunit'))
            wunit = xunit
        assert medium in ['air', 'vacuum', 'default']
    
        # Plot lines in spectrum medium (air/vacuum) if available
        if 'medium' in self.conditions and medium == 'default':
            medium = self.conditions['medium']

        def get_overlay(overlay):
            ''' Overlay line survey with a spectral quantity (like radiance or 
            transmittance) '''

            if isinstance(overlay, string_types):  # either get it from the Spectrum
                if overlay not in self.get_vars():
                    raise AttributeError('{0} not in variables list: {1}'.format(overlay,
                                         self.get_vars()))
                w, I = self.get(overlay, wunit=wunit, medium=medium)
                name = overlay
                units = self.units[overlay]
                return (w, I, name, units)

            else:  # Or use a given tuple or arrays
                try:
                    (w,I) = overlay
                except:
                    raise ValueError('Overlay has to be string, or (w,I) tuple of '+\
                                     'arrays')
                return (w, I, '', '')
            
            
        if overlay is not None:
            
            if type(overlay) is not list:
                overlay = [overlay]
            
            overlay = [get_overlay(ov) for ov in overlay]

            return LineSurvey(self, overlay=overlay, wunit=wunit, 
                              cutoff=cutoff, medium=medium, 
                              *args, **kwargs)

        else:
            return LineSurvey(self, wunit=wunit, cutoff=cutoff, medium=medium, 
                              *args, **kwargs)



    def get_conditions(self):
        ''' Get all physical / computational parameters.
        '''

        return self.conditions


    def print_conditions(self):
        ''' Prints all physical / computational parameters.
        
        '''
        
        return print_conditions(self.get_conditions(), self.cond_units)

    def store(self, path, discard=['lines', 'populations'], compress=False, 
              add_info=None, add_date=None, if_exists_then='error', verbose=True):
        ''' Save a Spectrum object in JSON format. Object can be recovered with
        :func:`~radis.tools.database.load_spec`. If many Spectrum are saved in a 
        same folder you can view their properties with the :class:`~radis.tools.database.SpecDatabase`
        structure.
    
        Parameters    
        ----------

        path: path to folder (database) or file
            if a folder, file is saved to database and name is generated automatically.
            if a file name, then Spectrum is saved to this file and the later
            formatting options dont apply

        file: str
            explicitely give a filename to save
    
        compress: boolean
            if True, removes all quantities that can be regenerated with the 
            :meth:`~radis.spectrum.spectrum.Spectrum.update` method
            e.g, transmittance if abscoeff and path length are given, radiance if
            emisscoeff and abscoeff are given in non-optically thin case, etc.
            Default False

        add_info: list
            append these parameters and their values if they are in conditions
            example::

                add_info = ['Tvib', 'Trot']

        discard: list of str
            parameters to exclude, for instance to save some memory for instance
            Default [`lines`, `populations`]: retrieved Spectrum looses the 
            :meth:`~radis.spectrum.spectrum.Spectrum.line_survey` ability, 
            and :meth:`~radis.spectrum.spectrum.Spectrum.plot_populations`
            (but it saves tons of memory!)

        if_exists_then: 'increment', 'replace', 'error'
            what to do if file already exists. If increment an incremental digit
            is added. If replace file is replaced (yeah). If error (or anything else)
            an error is raised. Default `error`


        Returns
        -------

        Returns filename used


        Notes
        -----
        
        If many spectra are stored in a folder, it may be time to set up a 
        :class:`~radis.tools.database.SpecDatabase` structure to easily see all 
        Spectrum conditions and get Spectrum that suits specific parameters. 
        

        Implementation:
            
            Shouldnt rely on a Database. One may just want to store/load a Spectrum
            once.
            
            
        Examples
        --------
        
        Store a spectrum in compressed mode, regenerate quantities after loading::
        
            from radis import load_spec 
            s.store('test.spec', compress=True)   # s is a Spectrum 
            s2 = load_spec('test.spec')
            s2.update()                           # regenerate missing quantities 
            
        
        See Also
        --------
        
        :class:`~radis.tools.database.SpecDatabase`, 
        :func:`~radis.tools.database.load_spec`, 
        :meth:`~radis.spectrum.spectrum.Spectrum.store`
            
        '''
        # Todo: maybe move most of the code here (filename writing) in database.py ?

        from radis.tools.database import save

        if isinstance(discard, string_types):
            discard = [discard]

        return save(self, path, discard=discard, compress=compress, add_info=add_info, 
                    add_date=add_date, if_exists_then=if_exists_then, verbose=verbose)

    def save(self, *args, **kwargs):
        ''' Alias to Spectrum.store. See Spectrum.store for documentation '''

        return self.store(*args, **kwargs)

    def resample(self, w_new, unit='same', out_of_bounds='nan', 
                 if_conflict_drop='error', medium='default', 
                 energy_threshold=1e-3, print_conservation=False):
        ''' Resample spectrum over a new wavelength
        Fills with transparent medium when out of bound (transmittance 1,
        radiance 0)

        Warning. This may result in information loss. Resampling is done with
        oversampling and spline interpolation. These parameters can be adjusted,
        and energy conservation ensured with the appropriate parameters.
        
        Uses the :func:`~radis.misc.signal.resample` function. 

    
        Parameters    
        ----------

        w_new: array
            new wavespace to resample the spectrum on. Must be inclosed in the
            current wavespace (we won't extrapolate)

        unit: 'same', 'nm', 'cm-1'
            unit of new wavespace. It 'same' it is assumed to be the current
            waveunit. Default 'same'. The spectrum waveunit is changed to this
            unit after resampling (i.e: a spectrum calculated and stored in `cm-1`
            but resampled in `nm` will be stored in `nm` from now on). 
            If 'nm', see also argument ``medium``
            
        out_of_bounds: 'transparent', 'nan', 'error'
            what to do if resampling is out of bounds. 'transparent': fills with
            transparent medium. 'nan': fill with nan. 'error': raises an error.
            Default 'nan'

        if_conflict_drop: 'error', 'convoluted', 'non_convoluted'
            There is a problem if both convoluted and non convoluted (*no_slit)
            quantities coexists, as they aren't scaled on the same wavespace
            grid. If 'error' an error is raised. If 'convoluted', convoluted 
            quantities will be dropped. If 'non_convoluted' non convoluted quantities 
            are dropped. Default 'error'

            TODO (but dangerous): reapply_slit at the end of the process if slit
            is in conditions?
            
        medium: 'air', 'vacuum', or 'default'
            in which medium is the new waverange is calculated if it is given 
            in 'nm'. Ignored if unit='cm-1'
            
            
        Other Parameters
        ----------------
        
        Inputs transfered to resample:

        energy_threshold: float
            if energy conservation (integrals) is above this threshold, raise an
            error.

        print_conservation: boolean
            if True, prints energy conservation. Default False.

            
        See Also
        --------
        
        :func:`~radis.misc.signal.resample`
            
        '''

        # Check inputs (check for deprecated)

        # ... see if convoluted / non convoluted values co-exist
        if 'wavespace' in self._q and 'wavespace' in self._q_conv:
            try:
                assert len(self._q['wavespace'])==len(self._q_conv['wavespace'])
                assert np.allclose(self._q['wavespace'], self._q_conv['wavespace'])
            except AssertionError:  # wavespaces are not the same.
                if if_conflict_drop == 'convoluted':
                    for q in list(self._q_conv.keys()):
                        del self._q_conv[q]
                elif if_conflict_drop == 'non_convoluted':
                    for q in list(self._q.keys()):
                        del self._q[q]
                elif if_conflict_drop == 'error':
                    raise ValueError('Cant resample as there are convoluted and non '+\
                                 'convoluted quantities in the Spectrum object (and '+\
                                 'wavespace are not the same). Use '+\
                                 "`if_conflict_drop='convoluted' or 'non_convoluted'`")
                else:
                    raise ValueError('Unknown value for if_conflict_drop: {0}'.format(
                            if_conflict_drop))

        # ... assert all values have the same wavespace
        # Now true by construction
#        quantities = list(self.values())
#        for q in quantities[1:]:
#            if not allclose(q[0], quantities[0][0]):
#                raise ValueError('Not all wavespaces are the same. Cant resample')

        # Get wavespace units
        waveunit = self.get_waveunit()   # spectrum unit
        if unit == 'same':               # resampled unit
            unit = waveunit
        else:
            unit = cast_waveunit(unit)
        
        # Get current waverange in output unit, output medium   -> w 
        if unit == 'nm':         # if asking for wavelength, check in which medium (air or vacuum?)
            # get defaults
            current_medium = self.get_medium()
            if medium == 'default':
                medium = current_medium 
            # convert if needed
            if medium == current_medium:
                pass # correct medium already, do nothing
            else:
                # update medium
                self.conditions['medium'] = medium
                
            w = self.get_wavelength(medium=medium)   # in correct medium
        elif unit == 'cm-1':
            if medium != 'default':
                raise ValueError('resampling to cm-1 but `medium` is given. It '+\
                                 'wont change wavenumbers, so just dont use it')
            w = self.get_wavenumber()              
        else:
            raise ValueError('Unknown unit: {0}'.format(unit))
            
        # Update waveunit to new unit 
        if unit != waveunit:
            self.conditions['waveunit'] = unit

        # Get wavespace
        update_q = 'wavespace' in self._q
        update_q_conv = 'wavespace' in self._q_conv

        # Now let's resample
        def get_filling(variable):
            ''' Get out of bounds values for spectral quantity `variable` '''
            if out_of_bounds == 'transparent':
                # Fill with optically transparent medium
                if variable in ['transmittance', 'transmittance_noslit']:
                    fill_with = 1
                else:
                    fill_with = 0
            elif out_of_bounds == 'nan':
                fill_with = 'nan'
            elif out_of_bounds == 'error':
                fill_with = 'error'
            else:
                raise ValueError('Unexpected value for out_of_bound: {0}'.format(
                        out_of_bounds))
            return fill_with
        
        
        # There are different cases depending on the unit of w_new
        # ... Note for devs: we're looping over dictionaries directly rather than
        # ... using the (safter) .get() function because it's much faster (the 
        # ... air2vacuum conversion in particular is quite slow, but has been 
        # ... done once for all with get_wavelength() above )
        if update_q:
            for (k, I) in self._q.items():
                if k == 'wavespace':
                    continue
                fill_with = get_filling(k)
                Inew = resample(w, I, w_new, ext=fill_with,
                                energy_threshold=energy_threshold,
                                print_conservation=False)
                self._q[k] = Inew
            # update wavespace
            self._q['wavespace'] = w_new

        if update_q_conv:
            for (k, I) in self._q_conv.items():
                if k == 'wavespace':
                    continue
                fill_with = get_filling(k)
                Inew = resample(w, I, w_new, ext=fill_with,
                                energy_threshold=energy_threshold,
                                print_conservation=False)
                self._q_conv[k] = Inew
            # update wavespace
            self._q_conv['wavespace'] = w_new


    # %% ======================================================================
    # Semi public functions
    # ----------------
    # Access object properties (but don't manipulate the spectrum itself)
    # XXX =====================================================================

    def wavespace(self):
        ''' Returns whether this spectrum is defined in wavelength (nm) or
        wavenumber (cm-1) '''

        warn(DeprecationWarning('wavespace replaced with get_waveunit'))

        return self.get_waveunit()

    def get_waveunit(self):
        ''' Returns whether this spectrum is defined in wavelength (nm) or
        wavenumber (cm-1) '''

        return self.conditions['waveunit']

    def is_at_equilibrium(self):
        ''' Returns whether this spectrum is at equilibrium. The following properties
        are checked:
            
            Tvib = Trot = Tgas
            self_absorption = True
            
        '''
        
        cond = self.conditions
        
        try:
            assert cond['Tvib'] == cond['Tgas']
            assert cond['Trot'] == cond['Tgas']
            if 'overpopulation' in cond:
                assert cond['overpopulation'] is None
            assert cond['self_absorption']  # is True
            
            return True
            
        except AssertionError:
            return False 
        except KeyError as err:
            warn('Condition missing to know if spectrum is at equilibrium: {0}'.format(err)+\
                 '. Assumed False')
            return False  
        
    def get_medium(self):
        ''' Returns in which medium the spectrum is calculated (air or vacuum), 
        based on the value on the self_absorption key in conditions. If not given, raises an error'''
        
        try:
            return self.conditions['medium']
        except KeyError:
            raise KeyError("We need to know if Spectrum is calculated in air or vacuum, but "+\
                           "`medium` is not defined in conditions. Please add the "+\
                             "value manually")
        
    def is_optically_thin(self):
        ''' Returns whether the spectrum is optically thin, based on the value
        on the self_absorption key in conditions. If not given, raises an error'''
        
        try:
            return not self.conditions['self_absorption']
        except KeyError:
            raise KeyError("We need to know if Spectrum is optically thin, but "+\
                           "`self_absorption` is not defined in conditions. Please add the "+\
                             "value manually")
        
    def copy(self, copy_lines=True):
        ''' Returns a copy of this Spectrum object (performs a smart deepcopy) '''
        try:
            return self.__copy__(copy_lines=copy_lines)
        except MemoryError:
            raise MemoryError("during copy of Spectrum. If you don't need them, "+\
                              "droping lines before copying may save a lot of space: "+\
                              "del s.lines ; or, use copy_lines=False")

    def __copy__(self, copy_lines=True):
        ''' Generate a new spectrum object

        Note: using deepcopy would work but then the Spectrum object would be pickled
        and unpickled again. It's a little faster here


        Notes
        -----
        
        Performance:

            deepcopy: 3.32 ms
            initial : 10 ms
            no check test, optimised: 1.8 ms
            ... asymptote: without evenly spaced check, without copies: 1.84 ms

        '''

        # quantities = {s:(v[0].copy(), v[1].copy()) for (s,v) in self.items()}  #╪ 1.8 ms
#        quantities = dict(self.items())   # 912 ns, not a copy but no need as
#                                        # Spectrum() recreates a copy anyway
        quantities = dict(self._get_items())

        try:
            units = deepcopy(self.units)
        except AttributeError:
            units = None

        conditions = deepcopy(self.conditions)  # 39 us
        try:
            cond_units = deepcopy(self.cond_units)
        except AttributeError:
            cond_units = None

        lines = None
        if copy_lines:
            try:
                lines = self.lines.copy(deep=True)   # 143 us  (2 ms with deepcopy(lines))
            except AttributeError:
                pass

        try:
            populations = self.populations
        except AttributeError:
            populations = None

        waveunit = self.get_waveunit()    # 163 ns
        name = self.name

        # Generate copied Spectrum
        s = Spectrum(      # 1.51 ms
                quantities=quantities,
                conditions=conditions,
                cond_units=cond_units,
                populations=populations,
                lines=lines,
                units=units,
                waveunit=waveunit,
                name=name,
                warnings=False,   # saves about 3.5 ms on the Performance test object
                )
        
        # Add slit information
        try:
            wslit, Islit = self.get_slit()
            s._slit['wavespace'] = wslit  # in 'waveunit'
            s._slit['intensity'] = Islit
        except KeyError:
            # no slit to apply
            pass

        return s

    def compare_with(self, other, spectra_only=False, plot=True, wunit='default',
                     verbose=True, rtol=1e-5, ignore_nan=False, ignore_outliers=False,
                     normalize=False, **kwargs):
        ''' Compare Spectrum with another spectrum

    
        Parameters    
        ----------
        
        other: type Spectrum
            another Spectrum to compare with

        spectra_only: boolean, or str
            if True, only compares spectral quantities (in the same waveunit)
            and not lines or conditions. If str, compare a particular quantity
            name. If False, compare everything (including lines and conditions
            and populations). Default False

        plot: boolean
            if True, use plot_diff to plot all quantities for the 2 spectra
            and the difference between them. Default False.
            
        wunit: 'nm', 'cm-1', 'default'
            in which wavespace to compare (and plot). If default, natural wavespace
            of first Spectrum is taken

        rtol: float
            relative difference to use for spectral quantities comparison

        ignore_nan: boolean
            if True, nans are ignored when comparing spectral quantities
            
        ignore_outliers: boolean, or float
            if not False, outliers are discarded. i.e, output is determined by::
            
                out = (~np.isclose(I, Ie, rtol=rtol, atol=0)).sum()/len(I) < ignore_outliers
        
        normalize: bool
            Normalize the spectra to be ploted 

        Other Parameters
        ----------------
        
        kwargs: dict
            arguments are forwarded to :func:`~radis.spectrum.compare.plot_diff`

        Returns
        -------
        
        equals: boolean
            return True if spectra are equal (respective to tolerance defined by 
            rtol and other input conditions)
            
            
        Examples
        --------
        
        Compare two Spectrum objects, or specifically the transmittance::
        
            s1.compare_with(s2)
            s1.compare_with(s2, 'transmittance')

        '''
        
        from radis.spectrum.compare import plot_diff
        
        # Check inputs
        if not 0<=ignore_outliers<1:
            raise ValueError('ignore_outliers should be < 1, or False')
        if not is_spectrum(other):
            raise TypeError('2nd object is not a Spectrum: got class {0}'.format(
                    other.__class__))
        if isinstance(spectra_only, string_types):   # case where we compare all quantities
            if not spectra_only in self.get_vars():
                raise ValueError('{0} is not a spectral quantity in our Spectrum ({1})'.format(
                        spectra_only, self.get_vars()))
            if not spectra_only in other.get_vars():
                raise ValueError('{0} is not a spectral quantity in the other Spectrum ({1})'.format(
                        spectra_only, other.get_vars()))
        if verbose:  # print conditions
            what = spectra_only if isinstance(spectra_only, string_types) else 'all quantities'
            msg = 'compare {0} with rtol={1}'.format(what, rtol)
            if ignore_nan: msg += ', ignore_nan'
            if ignore_outliers: msg += ', ignore_outliers={0}'.format(ignore_outliers)
            print(msg)
        if not plot and len(kwargs)>0:
            raise ValueError('Unexpected argument: {0}'.format(kwargs))

        if wunit == 'default':
            wunit = self.get_waveunit()

        def _compare_dataframes(df1, df2, name):
            ''' 
    
            Parameters    
            ----------
            
            df1, df2: pandas Dataframe
                lines, or vib/rovib levels dataframes 
            
            name: str
                for error message
            '''

#            if compare_lists(df1.keys(), df2.keys(), verbose=False) != 1:
#                if verbose: print('... keys in {0} dont match:'.format(name))
#                compare_lists(list(df1.keys()), list(df2.keys()), 
#                              verbose=True)
#                out = False
#            elif compare_lists(df1.index, df2.index, verbose=False) != 1:
#                if verbose: print('... index in {0} dont match:'.format(name))
#                compare_lists(list(df1.index), list(df2.index), 
#                              verbose=True)
#                out = False
#            else:
#                out = (df1 == df2).all().all()
#                
#            return out
            
            from pandas.util.testing import assert_frame_equal
            
            try:
                assert_frame_equal(df1.sort_index(axis=0).sort_index(axis=1), 
                                   df2.sort_index(axis=0).sort_index(axis=1),
                                   check_names=True)
                out = True
            
            except AssertionError as err:
                if verbose: 
                    print('Comparing ', name)
                    print(err.args[0])
                out = False
                
            return out

        def _compare_variables(I, Ie):
            ''' Compare spectral quantities I and Ie '''
                        
            if ignore_nan:
                b = ~(np.isnan(I) + np.isnan(Ie))
                I = I[b]
                Ie = Ie[b]
            
            if ignore_outliers:
                out = (~np.isclose(I, Ie, rtol=rtol, atol=0)).sum()/len(I) < ignore_outliers
            else:
                out = np.allclose(I, Ie, rtol=rtol, atol=0)
                
            return bool(out)

        b = True
        if isinstance(spectra_only, string_types): # compare this quantity
            vars = [spectra_only]
        else:                    # compare all quantities
            b = set(self.get_vars()) == set(other.get_vars())
            if not b and verbose: print('... list of quantities dont match: {0} vs {1}'.format(
                    self.get_vars(), other.get_vars()))
            vars = [k for k in self.get_vars() if k in other.get_vars()]

        if spectra_only:
            # Compare spectral variables
            # -----------
            for k in vars:
                w, q = self.get(k, wunit=wunit)
                w0, q0 = other.get(k, wunit=wunit)
                if len(w) != len(w0):
                    print('Wavespaces have different length (in {0})'.format(k))
                    b1 = False
                else:
                    b1 = allclose(w, w0, rtol=rtol, atol=0)
                    b1 *= _compare_variables(q, q0)
                    if not b1 and verbose:
                        error = np.nanmax(abs(q/q0-1))
                        avgerr = np.nanmean(abs(q/q0-1))
                        print('...', k, 'dont match (up to {0:.1f}% diff.,'.format(
                                error*100)+' average {0:.1f}%)'.format(avgerr*100))
                b *= b1
                
                if plot:
                    try:
                        plot_diff(self, other, var=k, wunit=wunit, 
                                  normalize=normalize, verbose=verbose,
                                  **kwargs)
                    except:
                        print('... couldnt plot {0}'.format(k))

        else:
            # Compare spectral variables
            # -----------
            for k in vars:
                w, q = self.get(k, wunit=wunit)
                w0, q0 = other.get(k, wunit=wunit)
                if len(w) != len(w0):
                    print('Wavespaces have different length (in {0})'.format(k))
                    b1 = False
                else:
                    b1 = allclose(w, w0, rtol=rtol, atol=0)
                    b1 *= _compare_variables(q, q0)
                    if not b1 and verbose:
                        error = np.nanmax(abs(q/q0-1))
                        print('...', k, 'dont match (up to {0:.1f}% diff.)'.format(
                                error*100))
                b *= b1

                if plot:
                    try:
                        plot_diff(self, other, var=k, wunit=wunit, 
                                  normalize=normalize, verbose=verbose,
                                  **kwargs)
                    except:
                        print('... couldnt plot {0}'.format(k))

            # Compare conditions and units
            # -----------
            b1 = self.conditions == other.conditions
            b2 = self.cond_units == other.cond_units
            b3 = self.units == other.units
            if not b1 and verbose: print('... conditions dont match')
            if not b2 and verbose: print('... conditions units dont match')
            if not b3 and verbose: print('... units dont match')
            b *= b1 * b2 * b3

            # Compare lines
            # -----------
            if self.lines is None and other.lines is None:
                b4 = True
            elif self.lines is None:
                b4 = False
            elif other.lines is None:
                b4 = False
            else:
                b4 = _compare_dataframes(self.lines, other.lines, 'lines')
            if not b4 and verbose: print('... lines dont match')
            b *= b4

            # Compare populations
            # -----------
            if self.populations is None and other.populations is None:
                b5 = True
            elif self.populations is None:
                b5 = False
            elif other.populations is None:
                b5 = False
            else:
                # Compare keys in populations 
                b5 = True
                if compare_lists(self.populations, 
                                other.populations, 
                                verbose='if_different') == 1:
                    # same molecules, compare isotopes
                    for molecule, isotopes in self.populations.items():
                        if compare_lists(isotopes, 
                                        other.populations[molecule], 
                                        verbose='if_different') == 1:
                            # same isotopes, compare electronic states
                            for isotope, states in isotopes.items():
                                if compare_lists(states, 
                                                other.populations[molecule][isotope], 
                                                verbose='if_different') == 1:
                                    # same electronic states, compare levels + other information
                                    for state, content in states.items():
                                        for k, v in content.items():
                                            if k in ['vib', 'rovib']:
                                                b5 *= _compare_dataframes(
                                                        v,
                                                        other.populations[molecule][isotope][state][k],
                                                        'populations of {0}({1})(iso{2})'.format(
                                                                molecule, state, isotope))
                                            else:
                                                b5 *= (v == other.populations[molecule][isotope][state][k])
                                else:
                                    b5 = False
                                    if verbose: print('{0}(iso{1}) states are different (see above)'.format(
                                                molecule, isotope))
                        else:
                            b5 = False
                            if verbose: print('{0} isotopes are different (see above)'.format(molecule))
                else:
                    b5 = False
                    if verbose: print('Molecules are different (see above)')
            if not b5 and verbose: print('... populations dont match (see detail above)')
            b *= b5


            # Compare slit
            # -----------
            
            if len(self._slit) == len(other._slit) == 0:
                # no slit anywhere
                b6 = True
            elif len(self._slit) != len(other._slit):
                b6 = False
                if verbose: print('A spectrum has slit function array but the other doesnt')
            else:
                ws, Is = self.get_slit()
                ws0, Is0 = other.get_slit()
                if len(ws) != len(ws0):
                    if verbose: print('Slits have different length')
                    b6 = False
                else:
                    b6 = allclose(ws, ws0, rtol=rtol, atol=0)
                    b6 *= _compare_variables(Is, Is0)
                    if not b6 and verbose:
                        print('Slit functions dont match')
            b *= b6

        return bool(b)

    # %% ======================================================================
    # Private functions
    # ----------------
    # XXX =====================================================================

    def _init_annotations(self):
        ''' Annotations are used to give typing hints for get() and plot() functions,
        based on what quantities are available in the Spectrum object '''

        from radis.tools.slit import SLIT_SHAPES

        try: # Python >3.6 only
            self.get.__annotations__['var']=[]
            self.plot.__annotations__['var']=[]
            self.apply_slit.__annotations__['shape'] = SLIT_SHAPES

        except AttributeError:
            pass  # old Python version


    def _add_quantity(self, name, w, I, warnings=True):
        ''' Add quantity. 
        
        Note: creates a copy of the input array
        '''

        assert(len(w)==len(I))

        def check_wavespace(w):
            ''' If warnings, check that array is evenly spaced. Returns a copy
            of input array.
            
            Note: this check takes a lot of time! 
            Is is not performed if warnings is False'''
            w = np.array(w)              # copy 
            if warnings:
                # Check Wavelength/wavenumber is evently spaced
                if not evenly_distributed(w, tolerance=1e-5):
                    warn('Wavespace is not evenly spaced ({0:.3f}%) for {1}.'.format(
                            np.abs(np.diff(w)).max()/w.mean()*100, name)\
                         + ' This may create problems when convolving with slit function')
            return w

        if name in CONVOLUTED_QUANTITIES:
            # Add wavespace
            if 'wavespace' in self._q_conv:
                if not np.allclose(w, self._q_conv['wavespace']):
                    raise ValueError('wavespace for {0} doesnt correspond to existing wavespace'.format
                                     (name)+' for convoluted quantities')
            else:
                self._q_conv['wavespace'] = check_wavespace(w)   # copy

            # Add quantity itself
            self._q_conv[name] = np.array(I)              # copy

        elif name in NON_CONVOLUTED_QUANTITIES:
            # Add wavespace
            if 'wavespace' in self._q:
                if not np.allclose(w, self._q['wavespace']):
                    raise ValueError('wavespace for {0} doesnt correspond to existing wavespace'.format
                                     (name)+' for non convoluted quantities')
            else:
                self._q['wavespace'] = check_wavespace(w)   # copy

            # Add quantity itself
            self._q[name] = np.array(I)              # copy

        else:
            raise ValueError('Unknown quantity: {0}'.format(name))

        # also make the quantity accessible with s.[name] like Pandas dataframes (removed eventually)
       # setattr(self, name, quantity)   # Warning this makes another copy of it (it's a tuple!)

        # add to annotations   (Python >3.6)
        try:
            self.get.__annotations__['var'].append(name)
            self.plot.__annotations__['var'].append(name)
        except AttributeError:  # old Python version
            pass

    def __eq__(self, other):
        """Override the default Equals behavior"""
        return self.compare_with(other, verbose=False, plot=False)

    def __ne__(self, other):
        """Define a non-equality test"""
        return not self.__eq__(other)

    def __dir__(self):
        ''' Names shown with tab completion: remove certain attributes to simplify
        the use of this class (@minou).
        '''

#        attrs = super(Spectrum, self).__dir__()
        attrs = dir(type(self))     # Python 2 and 3 compatible
        exclude = ['clear', 'fromkeys', 'items', 'pop', 'popitem', 'setdefault',
                   'values']

        return [k for k in attrs if not k in exclude]
    
    def __str__(self):
        ''' Print all Spectrum attributes'''
        
        # Print name
        print('Spectrum Name: ', self.get_name())
        
        # Print spectral quantities
        print('Spectral Quantities')
        print('-'*40)
        for k, v in self._get_items().items():
            print(' '*2, k, '\t({0} points)'.format(len(v[0])))
            
        # Print populations
        print('Populations Stored')
        print('-'*40)
        try:
            for k, v in self.populations.items():
                print(' '*2, k, '\t\t', list(v.keys()))
        except:
            pass
        
        # Print conditions
        self.print_conditions()
        
        return '' #self.print_conditions()
    
    # the following is so that json_tricks.dumps and .loads can be used directly,
    # ie.:  
    # >>> s == loads(s.dumps())
    # Does not work yet.
    # see https://github.com/mverleg/pyjson_tricks#class-instances
#    def __json_encode__(self):
#        ''' Called by json_tricks.dumps '''
#        from radis.tools.database import _format_to_jsondict
#        return _format_to_jsondict(self, discard=[], compress=[], verbose=True)
#    
#    def __json_decode__(self, **attrs):
#        from radis.tools.database import _json_to_spec
#        print(attrs)
#        raise
#        return _json_to_spec(attrs)

# %% ======================================================================
# Test class function
# -------------------
# XXX =====================================================================

# Test class
def is_spectrum(a):
    ''' Returns whether a is a Spectrum object
    
    Parameters
    ----------
    
    a: anything
        a Python object
        
    Returns
    -------
    
    True if a is a Spectrum object
    
    Notes
    -----
    
    is_spectrum compares the object class name (str): in some cases the Spectrum 
    class gets imported twice (when databases are involved, mostly), and a purely 
    isinstance() comparison fails
    
    '''

    return (isinstance(a, Spectrum) or
             repr(a.__class__) == repr(Spectrum))


# Test functions
if __name__ == '__main__':
    from radis.test.spectrum.test_spectrum import _run_testcases
    print('Test spectrum: ', _run_testcases(debug=False))