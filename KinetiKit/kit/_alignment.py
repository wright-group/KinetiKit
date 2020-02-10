"""
Functions specializing on array shifting and alignment
"""

import sys
import os
import csv

import numpy as np
import matplotlib.pyplot as plt

def find_nearest(array, value):
    """
    Returns element and index of array that is nearest to value.
    
    Parameters
    ----------
    array : arr
        The array from which a value and index will be returned
    value : float
        The value that the returned element of `arr` will be nearest to
        
    Returns
    -------
    idx : int
        The index of `arr` corresponding to the nearest value to `value`
    nearest_val : float
        The nearest value of the array to `value`
    """
    
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    nearest_val = array[idx]
    return idx, nearest_val

def roll_by_array_shift(arrays, refarray, shift, direction='forward'):
    """
    Performs ``numpy.roll(arrays)`` by the same number of points that 
    correspond to a `shift` in a reference array (`refarray`).
    
    Assumes linear step size for `refarray` and `arrays`.
    
    Parameters
    ----------
    arrays : 1D or 2D array
        Arrays whose values will be shifted.
    refarray : 1D array
        Array to whose values the `shift` pertains. Must have a 1:1 element
        correspondence with rollarray, even though it can be smaller. 
        `rollarray` will be rolled by the same number of elements that add up 
        to `shift` in `refarray`
    shift : float
        Desired shift in `refarray` values.
    direction : str, 'forward' or 'back'
        Determines the direction of rolling. Default is 'forward.'
    
    Returns
    -------
    rolled_array : arr
        Shifted array
    """
    rolled_arrays = arrays.copy()
    
    offset = shift + refarray[0]
    idx, val = find_nearest(refarray, offset)
    if direction == 'forward':
        rolled_arrays = np.roll(rolled_arrays, idx, axis=-1)
    elif direction == 'back':
        rolled_arrays = np.roll(rolled_arrays, -idx, axis=-1)    
    else:
        raise ValueError('Direction must be \'forward\' or \'back\'.')
        #issue_error
    
    return rolled_arrays

 
def align_by_max(arrays, refarray_x, refarray_y=None, avgnum=1, value=0):
    """
    Determines how many indices refarray_y must be rolled by so that its max.
    region aligns with `value` in `refarray_x', and then rolls each array in 
    `arrays` by the same amount. The number of points considered to define the
    maximum "region" of refarray_y is defined by `avgnum`.
    
    Parameters
    ----------
    arrays : 1-D or 2-D array
    refarray_x : 1-D array
    
    Optional Parameters
    -------------------
    refarray_y : 1-D array
        Default is None, in which case refarray_y is defined as the first array
        in  `arrays`.
    avgnum : integer
        The integer `avgnum` is used as follows: The top `avgnum` elements of 
        each array_list are considered, which each correspond to a unique value
        in `refarray`. For each array to be rolled, the average of the 
        `refarray` values corresponding to the `avgnum` maxima is considered, 
        and ultimately determines the horizontal shift. Default is 1. Setting 
        `avgnum` > 1 can be useful if the y array is noisy and a single maximum
        element is not representative of its true maximum.
    value : float
        Default is zero. If the exact value does not exist in `refarray_x`, 
        then the nearest value on the array will be chosen.
    
    Returns
    -------
    rolled_arrays : list of arrays or 2D array
        Set of aligned arrays.
    """
    rolled_arrays = arrays.copy()
    
    if refarray_y is None:
        if rolled_arrays.ndim == 2:
            refarray_y = rolled_arrays[0].copy()
        elif rolled_arrays.ndim == 1:
            refarray_y = rolled_arrays.copy()
        else: 
            print('Function only accepts 1D and 2D arrays.')
            #issue_error
            return
    
    sacrificial = refarray_y.copy()
    refvalues_of_maxes = np.array([])
        
    for i in range(0,avgnum):
        current_max = sacrificial.max()
        max_idx = np.where(refarray_y==current_max)[0][0]
        current_refvalue = refarray_x[max_idx]
        refvalues_of_maxes = np.append(refvalues_of_maxes, [current_refvalue])
        sacrificial[np.where(sacrificial==current_max)[0][0]] = -sys.maxsize
        
    avg_ref_value = np.average(refvalues_of_maxes)
    max_idx = find_nearest(refarray_x, avg_ref_value)[0]
        
    val_idx, val = find_nearest(refarray_x, value)
        
    shift_idx = int(max_idx - val_idx)
    
    rolled_arrays = np.roll(rolled_arrays, -shift_idx, axis=-1)
    
    return rolled_arrays


def align_by_steep(arrays, refarray_x, refarray_y=None, avgnum=3, value=0):
    """
    Analogous function to  ``align_by_max``: Rolls each array in `arrays` by 
    the index required to align the steepest region of `refarray_y` with a 
    `value` in `refarray_x`. The number of points considered to define the
    steepest "region" of refarray_y is defined by `avgnum`.
    
    Required Parameters
    -------------------
    arrays : array or list of 1D arrays
    refarray_x : 1-D array
    
    Optional Parameters
    -------------------
    refarray_y : 1-D array
        Default is None, in which case refarray_y is defined as the first array
        in  `arrays`.
    value : float
        Default is zero. If the exact value does not exist in `refarray_x`, 
        then the nearest value on the array will be chosen.
    
    Returns
    -------
    rolled_arrays : list of arrays or 2D array
        Set of aligned arrays.
    """
     
    rolled_arrays = arrays.copy()
    
    if refarray_y is None:
        if rolled_arrays.ndim == 2:
            refarray_y = rolled_arrays[0].copy()
        elif rolled_arrays.ndim == 1:
            refarray_y = rolled_arrays.copy()
        else: 
            print('Function only accepts 1D and 2D arrays.')
            #issue_error
            return
        
    difs = np.zeros(len(refarray_y))
    for i in range(len(refarray_y)-1):
        difs[i] = refarray_y[i]-refarray_y[i-2*avgnum]
    difs = np.roll(difs, -avgnum)
    
    steep_idx = np.argmax(difs)
    val_idx, val = find_nearest(refarray_x, value)
    shift_idx = int(steep_idx - val_idx)
    
    rolled_arrays = np.roll(rolled_arrays, -shift_idx, axis=-1)
        
    return rolled_arrays