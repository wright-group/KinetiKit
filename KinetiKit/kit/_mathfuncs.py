"""
Simple mathematical methods used throughout the module
"""

import sys
import os
import csv

import numpy as np
import matplotlib.pyplot as plt

def normalized(ar, alert=True):
    """
    Returns an array normalized by its maximum value.
    
    Assumes that ar is to be normalized by a positive value.
    
    Parameters
    ----------
    ar : arr
        The array to be normalized
    
    alert : boolean, optional
        Determines whether or not to print a warning if the maximum of 
        ar_max is less than or equal to zero. Default is True.
    
    Returns
    -------
    norm_ar : arr
        Normalized array
    """
    norm_ar = ar.copy()
    
    if ar.ndim > 2:
        print('Error: Normalized function accepts only 2D or 1D arrays')
        #issue_error
        return
    
    if ar.ndim == 1:
        norm_ar = np.reshape(norm_ar, (1, len(ar)))
        
    for i, raw_ar in enumerate(norm_ar):
        raw_ar_max = raw_ar.max()
        if raw_ar_max <= 0:
            norm_ar[i] = raw_ar
            if alert:
                print('Array does not contain non-zero positive values!')
        else:  
            norm_ar[i] = raw_ar / raw_ar_max
            
    if ar.ndim == 1:
        return norm_ar[0]
    else:
        return norm_ar

def precision(number):
    """
    Returns the last decimal position of the number

    Parameters
    ----------
    number : float
        Number whose precision to determine. Note: if an integer with a 
        floating point (e.g `1.`) is input, it will be considered a 
        float number with precision 1.

    Returns
    -------
    precision : int
        The negative log ofo a float containing the digit 1 in the same 
        position as the last decimal point of the input number. If the number
        does not contain decimal points, then the value 0 is returned.
    
    
    """
    
    if '.' not in str(number):
        return 0
    

    decimals = str(number).split('.')[1]
    i=0
    while i < len(decimals):
        i += 1
        
    return i


def Gauss(t, a, t0, fwhm):
    """
    Returns a Gaussian shape given an array of x value and Gaussian parameters.
    
    Parameters
    ----------
    t : arr
        Array of x values
    a : float
        Amplitude of Gaussian
    t0 : float
        center of Gaussian shape on `t`
    fwhm : float
        Full-width at half-maximum of Gaussian shape
    
    Returns
    -------
    gauss : arr
        Array with same shape as t that follows a Gaussian curve
    """
    
    sigma = fwhm/(2*np.sqrt(2*np.log(2)))
    gauss = a * np.exp(-(t - t0)**2 / (2 * sigma**2))
    return gauss