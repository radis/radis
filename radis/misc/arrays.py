# -*- coding: utf-8 -*-
"""
Description
------------

Functions to deal with numpy arrays 

"""

from __future__ import absolute_import, division, print_function, unicode_literals


import numpy as np
from numpy import hstack
from scipy import interpolate
from scipy.interpolate import interp1d
from six.moves import map



# Normalize

def norm(a, normby=None, how='max'):
    ''' Normalize a numpy array with its maximum. Or normalize it with another
    vector. Works if array contains nans.
    
    
    Parameters    
    ----------
    
    normby: array, or None
        if array, norm with this other array's maximum. If None, normalize with
        its maximum.
    '''
    if normby is not None:
        normby = np.abs(normby)
    else:
        normby = np.abs(a)
    
    if how=='max':
        return a / np.nanmax(normby)
    elif how=='mean':
        return a / np.nanmean(normby)
    else:
        raise ValueError('Unknown normalisation method')

def norm_on(a, w, wmin=None, wmax=None, how='max'):
    ''' Normalize `a` on a specific range of `w` '''
    imin = np.argmin((w-wmin)**2) if wmin else None
    imax = np.argmin((w-wmax)**2) if wmax else None
    if imin is not None and imax is not None:
        if imin > imax:
            imin, imax = imax, imin
        elif imin == imax:
            imax += 1
    return norm(a, a[imin:imax], how=how)

def scale_to(a, b, k=1):
    ''' Scale function a to k*b '''
    return a * k * max(np.abs(b)) / max(np.abs(a))

def array_allclose(a, b, rtol=1e-5, atol=1e-8, equal_nan=True):
    ''' Returns wheter a and b are all close (element wise). If not the same size,
    returns False (instead of crashing like the numpy version). Cf numpy.allclose
    docs for more information. 
    
    
    Parameters    
    ----------
    
    a, b: arrays
    
    rtol: float
    
    atol: float
    
    equal_nan: bool
        whether to consider Nan's as equal. Contrary to the numpy version this
        one is set to True by default 
    
    '''
    
    if len(a) != len(b):
        return False
    
    return np.allclose(a, b, rtol=rtol, atol=atol, equal_nan=equal_nan)
    
def nantrapz(I, w, dx=1.0, axis=-1):
    ''' Returns np.nan(I, w) discarding nan '''
    b = ~np.isnan(I)
    return np.trapz(I[b], w[b], dx=dx, axis=axis)
    
# %%
#==============================================================================
# Numpy Function
#==============================================================================

def shift_array(t0, y0, shift, nmax=10, tab=0):
    '''     
    shift the array and interpolate when the shift is smaller than the time step

    
    Parameters    
    ----------
    
    t0: array-like
        x
    
    y0: array-like
        y
    
    shift: float
        value to shift x and y    
    
    nmax : int
        maximum interpolation step
    
    tab: int
        unless said otherwise (tab != 0), replace values with 0. Default 0 


    Returns
    -------
    
    t, y: array-like
        shifted arrays


    Note
    ----
    
    Only tested with constant timesteps    
    '''

    if shift == 0:
        return t0, y0

    t = np.copy(t0)
    y = np.copy(y0)

    # Get the interpolation step required
    dt = t[1] - t[0]
    n = dt / (abs(shift) % dt)
    n = round(n)

    if n > nmax:
        print('Interpolation step required', n, 'has been reduced to', nmax)
        n = nmax

    if n > 1:
        # Interpolate
        step = (max(t) - min(t)) / (len(t) - 1) / n
        tint = np.arange(min(t), max(t), step)
        tck = interpolate.splrep(t, y, s=0)
        y = interpolate.splev(tint, tck, der=0)
        t = tint

    # Find shift coordinates (don't forget the negative shift case)
    if shift > 0:
        try:
            n = np.argmax(t >= (t[0] + shift))
        except ValueError:  # Empty sequence (shift probably too important)
            n = len(t)
    else:
        try:
            n = np.argmax(t >= (t[-1] + shift))
        # Empty sequence (negative shift probably too important)
        except ValueError:
            n = 0

    # Shift
    if shift > 0:
        y[n:] = y[:-n]
        y[:n] *= 0
        y[:n] += tab
    else:
        y[:n] = y[-n:]
        y[n:] *= 0
        y[n:] += tab

    return t, y


def calc_diff(t1, v1, t2, v2):
    ''' Substract two vectors that may have slightly offset abscisses 
    interpolating the correct values 
    
    
    Parameters    
    ----------
    
    t1, v1: array_like
        first vector and its abscisses
        
    t2, v2: array_like
        second vector and its abscisses
        

    Returns
    -------
    
    tdiff, vdiff: array_like
        substracted vector and its abscisses
        
    '''
    
    t1, v1, t2, v2 = list(map(np.array, (t1, v1, t2, v2)))
    
    # Deal with inversely sorted case
    if t1[-1] < t1[0]:
        t1, v1 = t1[::-1], v1[::-1]
    if t2[-1] < t2[0]:
        t2, v2 = t2[::-1], v2[::-1]

    # Get the overlapping range
    b = np.logical_and(t2 > t1[0], t2 < t1[-1])

    tdiff = t2[b]
    v2 = v2[b]

    # Interpolate the correct values
    f = interp1d(t1, v1)
    v1 = f(tdiff)

    # Finally, substract:
    vdiff = v1 - v2

    return tdiff, vdiff


def find_nearest(array, searched, returnarg=False):
    ''' Return the closest elements in array of 'searched' array,Also returns a boolean index'''

    b = np.zeros_like(array, dtype=bool)

#    def find_nearest(array,value):
#    '''  assuming array is sorted. '''
#        idx = np.searchsorted(array, value, side="left")
#        print('idx',idx)
#        if math.fabs(value - array[idx-1]) < math.fabs(value - array[idx]):
#            return idx-1
#        else:
#            return idx

    def find_nearest(array, value):
        return (np.abs(array - value)).argmin()

    try:
        for s in searched:
            b[find_nearest(array, s)] = True
    except:
        b[find_nearest(array, searched)] = True

    if returnarg:
        out = array[b], b, np.argmax(b)
    else:
        out = array[b], b

    return out


def find_first(arr, treshold):
    ''' Return the index of the first element of the array arr whose value
    is more than the treshold '''

    return np.argmax(arr > treshold)


def autoturn(data, key=-1):
    ''' key value : 
        0 don't transpose
        1 : transpose
        -1 : auto : make sure the vectors are along the longest dimension
    '''

    if key == 0:
        return data
    elif key == 1:
        return data.transpose()
    elif key == -1:
        if data.shape[0] == max(data.shape):
            # Probably columns :
            return data.transpose()
        else:
            return data

def centered_diff(w):
    ''' Return w[i+1]-w[i-1]/2, same size as w'''
    dw = np.diff(w)
    return (hstack((dw, dw[-1])) + hstack((dw[0], dw)))/2
           
def evenly_distributed(w, tolerance=1e-5):
    ''' Make sure array `w` is evenly distributed
    
    
    Parameters    
    ----------
    
    w : numpy array
        array to test 
    
    tolerance: float
        absolute tolerance
    '''
    mean_step = np.diff(w).mean()
    return (np.abs((np.diff(w)-mean_step)) > tolerance).sum() == 0

def bining(I, ymin=None, ymax=None, axis=1):
    ''' Averages a I multi-dimensional array (typically an image) along the y axis
    bining(I) corresponds to I.mean(axis=1)
    Nan are not taken into account
    
    
    Parameters    
    ----------
    
    I: numpy array
        intensity
    
    ymin: int [0-I.shape[1]]
        If None, 0 is used. Default None.
    
    ymax: int [0-I.shape[1]]
        If None, I.shape[1] is used. Default None.
    
    axis: int
        Default 1
    '''
    I = np.array(I)   # convert to array in case it's a Pandas dataframe for instance
    if ymin is None:
        ymin = 0
    if ymax is None:
        ymax = I.shape[axis]
    if ymin < 0:
        print('Warning in bining. ymin ({0}) < 0'.format(ymin))
    if ymax > I.shape[axis]:
        print('Warning in bining. ymax ({0}) > yaxis length ({1})'.format(ymax, 
              I.shape[axis]))
    return np.nanmean(I[:, ymin:ymax], axis=axis)


def count_nans(a):
    ''' Nan are good but only in India '''
    
    return np.isnan(a).sum()


def logspace(xmin, xmax, npoints):
    ''' Returns points from xmin to xmax regularly distributed on a logarithm
    space. 
    Numpy's logspace does the same from 10**xmin to 10**xmax
    '''
    
    return np.logspace(np.log10(xmin), np.log10(xmax), npoints)