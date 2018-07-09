# -*- coding: utf-8 -*-
"""
Summary
-------

Functions to compare and plot comparison of :class:`~radis.spectrum.spectrum.Spectrum` 
objects


Routine Listings
----------------

- :func:`~radis.spectrum.compare.get_diff`
- :func:`~radis.spectrum.compare.get_distance`
- :func:`~radis.spectrum.compare.get_ratio`
- :func:`~radis.spectrum.compare.get_residual`
- :func:`~radis.spectrum.compare.get_residual_integral`
- :func:`~radis.spectrum.compare.plot_diff`


-------------------------------------------------------------------------------



"""

from __future__ import print_function, absolute_import, division, unicode_literals
from radis.misc.arrays import array_allclose
from radis.misc.curve import curve_substract, curve_distance, curve_divide
from radis.spectrum.spectrum import make_up, cast_waveunit
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.widgets import MultiCursor
import numpy as np
from publib import set_style, fix_style
from warnings import warn

# %% ======================================================================
# External functions
# ----------------
# XXX =====================================================================


def get_diff(s1, s2, var, wunit='default', Iunit='default', medium='default',
             resample=True):
    ''' Get the difference between 2 spectra
    Basically returns w1, I1 - I2 where (w1, I1) and (w2, I2) are the values of
    s1 and s2 for variable var. (w2, I2) is linearly interpolated if needed.

        .. math::

            dI = I_1 - I_2


    Parameters    
    ----------

    s1, s2: Spectrum objects

    var: str
        spectral quantity (ex: 'radiance', 'transmittance'...)

    wunit: 'nm', 'cm-1'
        waveunit to compare in

    Iunit: str
        if ``'default'`` use s1 unit for variable var

    medium: 'air', 'vacuum', default'
        propagating medium to compare in (if in wavelength)

    Returns    
    -------

    w1, Idiff: array
        difference interpolated on the second range(order?)
    
    Notes
    -----
    
    Uses :func:`~radis.misc.curve.curve_substract` internally

    See Also
    --------

    :func:`~radis.spectrum.compare.get_ratio`, 
    :func:`~radis.spectrum.compare.get_distance`,  
    :func:`~radis.spectrum.compare.get_residual`,
    :func:`~radis.spectrum.compare.get_residual_integral`, 
    :func:`~radis.spectrum.compare.plot_diff` 
    :meth:`~radis.spectrum.spectrum.compare_with` 
    '''

    w1, I1, w2, I2 = _get_defaults(s1, s2, var=var, wunit=wunit, Iunit=Iunit,
                                   medium=medium, assert_same_wavelength=not resample)

    # basically w1, I1 - I2 (on same range)
    return curve_substract(w1, I1, w2, I2)


def get_ratio(s1, s2, var, wunit='default', Iunit='default', medium='default',
              resample=True):
    ''' Get the ratio between two spectra
    Basically returns w1, I1 / I2 where (w1, I1) and (w2, I2) are the values of
    s1 and s2 for variable var. (w2, I2) is linearly interpolated if needed.

        .. math::

            R = I_1 / I_2


    Parameters    
    ----------

    s1, s2: Spectrum objects

    var: str
        spectral quantity 

    wunit: 'nm', 'cm-1'
        waveunit to compare in

    Iunit: str
        if ``'default'`` use s1 unit for variable var

    medium: 'air', 'vacuum', default'
        propagating medium to compare in (if in wavelength)

    Notes
    -----
    
    Uses :func:`~radis.misc.curve.curve_divide` internally

    See Also
    --------

    :func:`~radis.spectrum.compare.get_diff`, 
    :func:`~radis.spectrum.compare.get_distance`,  
    :func:`~radis.spectrum.compare.get_residual`,
    :func:`~radis.spectrum.compare.get_residual_integral`, 
    :func:`~radis.spectrum.compare.plot_diff`,
    :meth:`~radis.spectrum.spectrum.compare_with` 
    

    '''

    w1, I1, w2, I2 = _get_defaults(s1, s2, var=var, wunit=wunit, Iunit=Iunit,
                                   medium=medium, assert_same_wavelength=not resample)

    # basically w1, I1 - I2 (on same range)
    return curve_divide(w1, I1, w2, I2)


def get_distance(s1, s2, var, wunit='default', Iunit='default', medium='default',
                 resample=True):
    ''' Get a regularized Euclidian distance between two spectra ``s1`` and ``s2`` 

    This regularized Euclidian distance minimizes the effect of a small shift in 
    wavelength between the two spectra

    .. math::

        D(w_1)[i] = \sqrt{ \sum_j (\hat{I_1}[i]  - \hat{I_2}[j] )^2 + (\hat{w_1}[i] - \hat{w_2}[j])^2}
        
    Where values are normalized as:
        
    .. math::

        \hat{A} = \\frac{A}{max(A) - min(A)}

    If waveranges dont match, ``s2`` is interpolated over ``s1``.

    .. warning :: 
        
        This is a distance on both the waverange and the intensity axis. 
        It may be used to compensate for a small offset in your experimental 
        spectrum (due to wavelength calibration, for instance) but can lead
        to wrong fits easily. Plus, it is very cost-intensive! Better use 
        :func:`~radis.spectrum.compare.get_residual` for an automatized procedure. 


    Parameters    
    ----------

    s1, s2: Spectrum objects

    var: str
        spectral quantity 

    wunit: 'nm', 'cm-1'
        waveunit to compare in

    Iunit: str
        if ``'default'`` use s1 unit for variable var

    medium: 'air', 'vacuum', default'
        propagating medium to compare in (if in wavelength)
    
    Notes
    -----
    
    Uses :func:`~radis.misc.curve.curve_distance` internally

    See Also
    --------

    :func:`~radis.spectrum.compare.get_diff`, 
    :func:`~radis.spectrum.compare.get_ratio`,  
    :func:`~radis.spectrum.compare.get_residual`,
    :func:`~radis.spectrum.compare.get_residual_integral`, 
    :func:`~radis.spectrum.compare.plot_diff`,
    :meth:`~radis.spectrum.spectrum.compare_with` 

    '''
    # TODO: normalize with Imax, wmax 

    w1, I1, w2, I2 = _get_defaults(s1, s2, var=var, wunit=wunit, Iunit=Iunit,
                                   medium=medium, assert_same_wavelength=not resample)

    # euclidian distance from w1, I1
    return curve_distance(w1, I1, w2, I2, discard_out_of_bounds=True)


def get_residual(s1, s2, var, norm='L2', ignore_nan=False):
    ''' Returns L2 norm of ``s1`` and ``s2``

    For ``I1``, ``I2``, the values of variable ``var`` in ``s1`` and ``s2``, 
    respectively, residual is calculated as:
        
    For ``L2`` norm:

        .. math::
    
            res = \\frac{\\sqrt{\\sum_i {(s_1[i]-s_2[i])^2}}}{N}.

    For ``L1`` norm:
            
        .. math::
    
            res = \\frac{\\sqrt{\\sum_i {|s_1[i]-s_2[i]|}}}{N}.


    Parameters    
    ----------

    s1, s2: :class:`~radis.spectrum.spectrum.Spectrum` objects
        if not on the same range, ``s2`` is resampled on ``s1``.

    var: str
        spectral quantity

    norm: 'L2', 'L1'
        which norm to use 

    Other Parameters
    ----------------

    ignore_nan: boolean
        if ``True``, ignore nan in the difference between s1 and s2 (ex: out of bound)
        when calculating residual. Default ``False``. Note: ``get_residual`` will still 
        fail if there are nan in initial Spectrum. 

    Notes
    -----

    0 values for I1 yield nans except if I2 = I1 = 0

    when s1 and s2 dont have the size wavespace range, they are automatically
    resampled through get_diff on 's1' range 

    Implementation of ``L2`` norm::

        np.sqrt((dI**2).sum())/len(dI)
    
    Implementation of ``L1`` norm::
        
        np.abs(dI).sum()/len(dI)


    See Also
    --------

    :func:`~radis.spectrum.compare.get_diff`, 
    :func:`~radis.spectrum.compare.get_ratio`, 
    :func:`~radis.spectrum.compare.get_distance`, 
    :func:`~radis.spectrum.compare.plot_diff`, 
    :func:`~radis.spectrum.compare.get_residual_integral`, 
    :meth:`~radis.spectrum.spectrum.compare_with` 
    '''

    w, I = s1.get(var)
    # mask for 0
    wdiff, dI = get_diff(s1, s2, var, resample=True)

    if ignore_nan:
        b = np.isnan(dI)
        wdiff, dI = wdiff[~b], dI[~b]

    if norm == 'L2':
        return np.sqrt((dI**2).sum())/len(dI)
    elif norm == 'L1':
        return (np.abs(dI)).sum()/len(dI)
    else:
        raise ValueError('unexpected value for norm')


def get_residual_integral(s1, s2, var, ignore_nan=False):
    ''' Returns integral of the difference between two spectra s1 and s2, 
    relatively to the integral of spectrum s1 
    
    Compared to :func:`~radis.spectrum.compare.get_residual`, this tends to 
    cancel the effect of the gaussian noise of an experimental spectrum. 
    
    ::
    
        res = trapz(I2-I1, w1) / trapz(I1, w1)

    Note: when the considered variable is ``transmittance`` or ``transmittance_noslit``, 
    the *upper* integral is used (up to 1) to normalize the integral difference
    
    ::
    
        res = trapz(I2-I1, w1) / trapz(1-I1, w1)

    Parameters
    ----------

    s1, s2: Spectrum objects

    var: str
        spectral quantity

    Other Parameters
    ----------------

    ignore_nan: boolean
        if ``True``, ignore nan in the difference between s1 and s2 (ex: out of bound)
        when calculating residual. Default ``False``. Note: ``get_residual_integral``
        will still fail if there are nan in initial Spectrum. 

    Notes
    -----

    For I1, I2, the values of 'var' in s1 and s2, respectively, residual
    is calculated as::

        res = trapz(I2-I1, w1) / trapz(I1, w1)

    0 values for I1 yield nans except if I2 = I1 = 0

    when s1 and s2 dont have the size wavespace range, they are automatically
    resampled through get_diff on 's1' range 


    See Also
    --------

    :func:`~radis.spectrum.compare.get_diff`, 
    :func:`~radis.spectrum.compare.get_ratio`, 
    :func:`~radis.spectrum.compare.get_distance`, 
    :func:`~radis.spectrum.compare.get_residual`, 
    :func:`~radis.spectrum.compare.plot_diff`,
    :meth:`~radis.spectrum.spectrum.compare_with` 
    '''

    w, I = s1.get(var)
    # mask for 0
    wdiff, dI = get_diff(s1, s2, var, resample=True)

    if ignore_nan:
        b = np.isnan(dI)
        wdiff, dI = wdiff[~b], dI[~b]
        b = np.isnan(I)
        w, I = w[~b], I[~b]
        
    if var in ['transmittance', 'transmittance_noslit']:
        norm = 1 - np.trapz(I, w)
    else:
        norm = np.trapz(I, w)

    return np.abs(np.trapz(dI, wdiff)) / norm

#    b1 = (I_avg == 0)
#    b2 = (dI == 0)

#    with catch_warnings():
#        filterwarnings('ignore', 'invalid value encountered in true_divide')
#        res = dI/I_avg

#    res[b1*b2] = 0             # if both I_avg and dI = 0, solve nan as 0

#    return np.sqrt((res**2).sum())/len(res)
#    return res.mean()


def _get_defaults(s1, s2, var, wunit='default', Iunit='default', medium='default',
                  assert_same_wavelength=False):
    ''' See get_distance, get_diff '''

    # Check inputs, get defaults
    # ----
    if Iunit == 'default':
        try:
            Iunit = s1.units[var]
        except KeyError:  # unit not defined in dictionary
            raise KeyError('Iunit not defined in spectrum. Cant plot')
    # Format units
    if wunit == 'default':
        wunit = s1.get_waveunit()
    wunit = cast_waveunit(wunit)
    if medium == 'default':
        medium = s1.conditions.get('medium', None)

    # Get data
    # ----
    w1, I1 = s1.get(var, wunit=wunit, Iunit=Iunit, medium=medium)
    w2, I2 = s2.get(var, wunit=wunit, Iunit=Iunit, medium=medium)

    if assert_same_wavelength:
        if not array_allclose(w1, w2):
            raise AssertionError(
                'Wavespace are not the same: use Spectrum.resample()')

    return w1, I1, w2, I2


def plot_diff(s1, s2, var=None, wunit='default', Iunit='default', medium='default',
              resample=True, method='diff',  show_points=False,
              label1=None, label2=None, figsize=None, title=None, nfig=None,
              normalize=False, verbose=True, save=False, show=True):
    ''' Plot two spectra, and the difference between them

    If waveranges dont match, ``s2`` is interpolated over ``s1``. 


    Parameters    
    ----------

    s1, s2: Spectrum objects

    var: str, or None
        spectral quantity to plot (ex: ``'abscoeff'``). If None, 
        plot the first one in the Spectrum from ``'radiance'``, 
        ``'radiance_noslit'``, ``'transmittance'``, etc.

    wunit: ``'default'``, ``'nm'``, ``'cm-1'``
        if ``'default'``, use first spectrum wunit

    Iunit: str
        if ``'default'``, use first spectrum unit

    medium: ``'air'``, ``'vacuum'``, ``'default'``
        if ``'default'``, use first spectrum propagating medium

    method: ``'distance'``, ``'diff'``, ``'ratio'``
        If ``'diff'``, plot difference at same wavespace position. 
        If ``'distance'``, plot Euclidian distance (note that units are meaningless then)
        If ``'ratio'``, plot ratio of two spectra
        Default ``'diff'``.

        .. warning::
            with ``'distance'``, calculation scales as ~N^2 with N the number
            of points in a spectrum (against ~N with ``'diff'``). This can quickly 
            override all memory.

    normalize: bool
        Normalize the spectra to be ploted 

    Other Parameters
    ----------------

    show_points: boolean
        if True, make all points appear with 'o'

    label1, label2: str
        curve names

    figsize
        figure size

    nfig: int, str
        figure number of name

    title: str
        title

    verbose: boolean
        if True, plot stuff such as rescale ratio in normalize mode. Default ``True``

    save: str
        Default is False. By default won't save anything, type the path of the 
        destination if you want to save it (format in the name).
        
    show: Bool
        Default is True. Will show the plots : bad if there are more than 20.
    Examples
    --------

    ::

        Punit = 'mW/cm2/sr'
        fig, axes = plot_diff(s10, s50, figsize=(18,6),
              label1='brd 10 cm-1, P={0:.2f} {1}'.format(s10.get_power(unit=Punit),Punit),
              label2='brd 50 cm-1, P={0:.2f} {1}'.format(s50.get_power(unit=Punit),Punit)
              )



    See Also
    --------

    :func:`~radis.spectrum.compare.get_diff`, 
    :func:`~radis.spectrum.compare.get_ratio`, 
    :func:`~radis.spectrum.compare.get_distance`, 
    :func:`~radis.spectrum.compare.get_residual`, 
    :meth:`~radis.spectrum.spectrum.compare_with` 

    '''

    # Get defaults
    # ---
    if var is None:    # if nothing is defined, try these first:
        params = s1.get_vars()
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
    if Iunit == 'default':
        try:
            Iunit = s1.units[var]
        except KeyError:  # unit not defined in dictionary
            raise KeyError('Iunit not defined in spectrum. Cant plot')
    if wunit == 'default':
        wunit = s1.get_waveunit()
    if medium == 'default':
        medium = s1.conditions.get('medium', None)

    # Get data
    # ----
    if normalize:
        w1, I1 = s1.get(var, wunit, Iunit, medium)
        w2, I2 = s2.get(var, wunit, Iunit, medium)
        if method == 'distance':
            wdiff, Idiff = curve_distance(w1, I1/np.max(I1), w2, I2/np.max(I2), discard_out_of_bounds=True)
        elif method == 'diff':
            wdiff, Idiff = curve_substract(w1, I1/np.max(I1), w2, I2/np.max(I2))
        elif method == 'ratio':
            wdiff, Idiff = curve_divide(w1, I1/np.max(I1), w2, I2/np.max(I2))
        else:
            raise ValueError('Unknown comparison method: {0}'.format(method))
    else:
        if method == 'distance':
            wdiff, Idiff = get_distance(
                s1, s2, var=var, wunit=wunit, Iunit=Iunit, medium=medium)
        elif method == 'diff':
            wdiff, Idiff = get_diff(
                s1, s2, var=var, wunit=wunit, Iunit=Iunit, medium=medium)
        elif method == 'ratio':
            wdiff, Idiff = get_ratio(
                s1, s2, var=var, wunit=wunit, Iunit=Iunit, medium=medium)
        else:
            raise ValueError('Unknown comparison method: {0}'.format(method))

    # Plot
    # ----

    # Format units
    wunit = cast_waveunit(wunit)
    if wunit == 'cm-1':
        xlabel = 'Wavenumber (cm-1)'
    elif wunit == 'nm':
        xlabel = 'Wavelength (nm)'

    # Init figure
    set_style('origin')
    fig = plt.figure(num=nfig, figsize=figsize)
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])
    ax1.get_shared_x_axes().join(ax0, ax1)

    # Plotting style
    if show_points:
        style = '-o'
    else:
        style = '-'

    # Get labels and names
    if label1 is None:
        label1 = s1.get_name()
    if label2 is None:
        label2 = s2.get_name()
    # Plot compared spectra
    if normalize:
        w1, I1 = s1.get(var, wunit, Iunit, medium)
        w2, I2 = s2.get(var, wunit, Iunit, medium)
        if verbose:
            print(('Rescale factor: '+str(np.max(I1)/np.max(I2))))
        ax0.plot(w1, I1/np.max(I1), ls=style, color='k', lw=3, label=label1)
        ax0.plot(w2, I2/np.max(I2), ls=style, color='r', lw=1, label=label2)
    else:
        ax0.plot(*s1.get(var, wunit, Iunit, medium),
                 ls=style, color='k', lw=3, label=label1)
        ax0.plot(*s2.get(var, wunit, Iunit, medium),
                 ls=style, color='r', lw=1, label=label2)

    Iunit = make_up(Iunit)  # cosmetic changes

    ax0.tick_params(labelbottom='off')
    if label1 is not None or label2 is not None:
        ax0.legend(loc='best')

    # Start to 0
    if var in ['radiance_noslit', 'radiance', 'abscoeff', 'absorbance']:
        ax0.set_ylim(bottom=0)

    # plot difference (sorted)
    b = np.argsort(wdiff)
    ax1.plot(wdiff[b], Idiff[b], style, color='k', lw=1)

    if method == 'diff':
        fig.text(0.09, 0.38, 'diff')
    elif method == 'distance':
        fig.text(0.09, 0.38, 'distance')
    elif method == 'ratio':
        fig.text(0.09, 0.38, 'ratio')

    # Write labels
    ax1.set_xlabel(make_up(xlabel))
    if normalize:
        fig.text(0.02, 0.5, 'Arb. Units',
                 va='center', rotation='vertical')
    else:
        fig.text(0.02, 0.5, ('{0} ({1})'.format(make_up(var), Iunit)),
                 va='center', rotation='vertical')

    # Set limits
    if method == 'diff':
        # symmetrize error scale:
        ymin, ymax = ax1.get_ylim()
        ymax = max(abs(ymin), abs(ymax))
        ax1.set_ylim(-ymax, ymax)
    elif method == 'distance':
        ax1.set_ylim(bottom=0)
    elif method == 'ratio':
        # auto-zoom on min, max, but discard first decile (case of spikes / divergences)
        Idiff_s = np.sort(Idiff)
        ax1.set_ylim(bottom=0, top=Idiff_s[len(
            Idiff_s)//10]+Idiff_s[-len(Idiff_s)//10])

    if title:
        fig.suptitle(title)

    # Fix format
    fix_style('origin', ax=ax0)
    fix_style('origin', ax=ax1)
    plt.tight_layout()
    if title:
        plt.subplots_adjust(left=0.15, top=0.92)
    else:
        plt.subplots_adjust(left=0.15)

    # Add cursors
    fig.cursors = MultiCursor(fig.canvas, (ax0, ax1),
                              color='r', lw=1, alpha=0.2, horizOn=False,
                              vertOn=True)
    if show:
        plt.show()
    if save:
        fig.savefig(save)

    return fig, [ax0, ax1]
''' Return the average distance between two spectra.
    It's important to note that if averageDistance(s1, s2)==0 then s1 = s2
    

   .. math::

      \\sqrt{\\sum {(u_i-v_i)^2 / V[x_i]}}.
'''


def averageDistance(s1, s2, var='radiance'):
    ''' Return the average distance between two spectra.
    It's important to note for fitting that if averageDistance(s1, s2)==0 
    then s1 = s2

    Parameters    
    ----------

    s1, s2: Spectrum objects
        spectra to be compared

    var: str, optional
        spectral quantity (ex: 'radiance', 'transmittance'...)

    Returns    
    -------

    distance: float
        Average distance as in the following equation:

        .. math::

            dist = \\frac{\\sqrt{\\sum_i {(s1_i-s2_i)^2}}}{N}.

    '''
    
    warn(FutureWarning('will be deprecated in future versions. Use get_residual(s1, s2, norm=L2) instead'))
    
    distArray = get_diff(s1, s2, var)
    distArray_Y2 = distArray[1]*distArray[1]
    N = np.size(distArray[1])
    distance = np.sqrt(np.sum(distArray_Y2))/N
    return distance


if __name__ == '__main__':

    from radis.test.spectrum.test_compare import _run_testcases
    print(('Test compare.py: ', _run_testcases(verbose=True)))
