# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 16:29:18 2018

@author: erwan

Operations on Curves, where a curve is a (w, I) array tuple

Similar to Origin Simple Curve Math operators

"""

from __future__ import absolute_import
import numpy as np
from scipy.interpolate import interp1d
from scipy.spatial.distance import cdist

def curve_distance(w1, I1, w2, I2, discard_out_of_bounds=True):
    ''' Get euclidian distance from curve (w1, I1) to curve (w2, I2)
    
    Euclidian distance minimizes the effect of a small shift in between the two
    curves in case of stiff curves (like a spectrum bandhead can be)
    
    No interpolation needed neither. 
    
    Distances for out of bounds values is set to nan
    
    
    Parameters    
    ----------
    
    w1, I1: array
        range and values for first curve
    
    w2, I2: array
        range and values for 2nd curve
    
    discard_out_of_bounds: boolean
        if True, distance for out of bound values is set to nan. Else, it 
        will be the distance from the last point.
    

    Returns
    -------
    
    w1, Idist: array
        minimal distance from I1 to I2, for each point in (w1, I1)
    '''
    
    dist = cdist(np.array((w1, I1)).T, np.array((w2, I2)).T).min(axis=1)
    
    # discard out of bound values
    if discard_out_of_bounds:
        b = np.logical_or(w1 < w2.min(), w1 > w2.max())
        dist[b] = np.nan
    
    return w1, dist

def curve_add(w1, I1, w2, I2, is_sorted=False, kind='linear'):
    ''' Add curve (w2, I2) from (w1, I1) 
    Linearly interpolates if the two ranges dont match. Fills out of bound 
    parameters with nan.
    
    Similar to Origin "Simple Curve Math Substract"
    
    
    Parameters    
    ----------
    
    w1, I1: array
        range and values for first curve
    
    w2, I2: array
        range and values for 2nd curve
    
    is_sorted: boolean
        (optional) if True, doesnt sort input arrays
        
    kind: str
        interpolation kind. Default 'linear'. See scipy.interpolate.interp1d
    

    Returns
    -------
    
    w1, Idiff: array
        difference interpolated on the second range
        
    '''
    
    I1, I2 = _curve_interpolate(w1, I1, w2, I2, is_sorted=is_sorted, kind=kind)
        
    return w1, I1 + I2
  
def curve_substract(w1, I1, w2, I2, is_sorted=False, kind='linear'):
    ''' Substracts curve (w2, I2) from (w1, I1) 
    Linearly interpolates if the two ranges dont match. Fills out of bound 
    parameters with nan.
    
    Similar to Origin "Simple Curve Math Substract"
    
    
    Parameters    
    ----------
    
    w1, I1: array
        range and values for first curve
    
    w2, I2: array
        range and values for 2nd curve
    
    is_sorted: boolean
        (optional) if True, doesnt sort input arrays
        
    kind: str
        interpolation kind. Default 'linear'. See scipy.interpolate.interp1d
    

    Returns
    -------
    
    w1, Idiff: array
        difference interpolated on the second range
        
    '''
    
    I1, I2 = _curve_interpolate(w1, I1, w2, I2, is_sorted=is_sorted, kind=kind)
        
    return w1, I1 - I2
  
def curve_multiply(w1, I1, w2, I2, is_sorted=False, kind='linear'):
    ''' Multiply curve (w2, I2) with (w1, I1) 
    Linearly interpolates if the two ranges dont match. Fills out of bound 
    parameters with nan.
    
    Similar to Origin "Simple Curve Math Substract"
    
    
    Parameters    
    ----------
    
    w1, I1: array
        range and values for first curve
    
    w2, I2: array
        range and values for 2nd curve
    
    is_sorted: boolean
        (optional) if True, doesnt sort input arrays
        
    kind: str
        interpolation kind. Default 'linear'. See scipy.interpolate.interp1d
    

    Returns
    -------
    
    w1, Idiff: array
        difference interpolated on the second range
        
    '''
    
    I1, I2 = _curve_interpolate(w1, I1, w2, I2, is_sorted=is_sorted, kind=kind)
        
    return w1, I1 * I2
  
def curve_divide(w1, I1, w2, I2, is_sorted=False, kind='linear'):
    ''' Divides curve (w1, I1) by (w2, I2)
    Linearly interpolates if the two ranges dont match. Fills out of bound 
    parameters with nan.
    
    Similar to Origin "Simple Curve Math Substract"
    
    
    Parameters    
    ----------
    
    w1, I1: array
        range and values for first curve
    
    w2, I2: array
        range and values for 2nd curve
    
    is_sorted: boolean
        (optional) if True, doesnt sort input arrays
        
    kind: str
        interpolation kind. Default 'linear'. See scipy.interpolate.interp1d
    

    Returns
    -------
    
    w1, Idiff: array
        difference interpolated on the second range
        
    '''
    
    I1, I2 = _curve_interpolate(w1, I1, w2, I2, is_sorted=is_sorted, kind=kind)
        
    return w1, I1 / I2
  
def _curve_interpolate(w1, I1, w2, I2, is_sorted=False, kind='linear'):
    
    # First sort both
    if not is_sorted:
        b = np.argsort(w1)
        w1, I1 = w1[b], I1[b]
        b = np.argsort(w2)
        w2, I2 = w2[b], I2[b]
    
    # resample if needed
    if not np.array_equal(w1, w2):
        f = interp1d(w2, I2, kind=kind, bounds_error=False)
        I2 = f(w1)
        
    return I1, I2
