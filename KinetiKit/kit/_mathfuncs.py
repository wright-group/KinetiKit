"""
Simple mathematical methods used throughout the module
"""

import sys
import os
import csv

import numpy as np
from scipy import special
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

def ExpGauss(t, t0, a, fwhm, tau, b=0):
    """ 
    Parameters
    ----------
    t : numpy array
        The time array on which to build the IRF
    t0 : float
        Mean of Gaussian profile (before exponential modification)
    a : float
        Amplitude of function 
    fwhm : float
        The FWHM that determines the variance of either type of curve. Must 
        have the same time units as `t`.
    tau : float
        Lifetime of decay of exponentially weighted Gaussian. Ignored if 
        `weighted` is set to False.
    b : float
        Offset of the final function, recommended as zero.
        
    """
    
    k = 1/tau # rate of exponential decay
    sigma = fwhm/(2*np.sqrt(2*np.log(2)))
    arg1 = k/2*(2*t0 + k*sigma**2 - 2*t)
    arg2 = (t0 + k*sigma**2 - t) / (np.sqrt(2)*sigma)
    
    one = np.exp(arg1)
    two = special.erfc(arg2)
    return a * k/2 * one * two + b

def GaussDiff(t, t0, a, fwhm, tau, b, c, weighted, tau_wt):
    """ 
    Creates n Exponentially weighted Gaussian IRF with a diffusion tail
    
    Parameters
    ----------
    t : numpy array
        The time array on which to build the IRF
    irf_type : string
        Type of IRF to be constructed; ```irf_type = 'Gauss'``` will return a 
        normal Gaussian curve, while ```irf_type = 'GaussDiff'``` will return
        an exponentially modified Gaussian with a diffusion tail.
    weighted : Boolean
        Whether to use a pure Gaussian or exponentially weighted Gaussian for 
        the main "peak" of the IRF. Default is True. 
        See en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution. 
    fwhm : float
        The FWHM that determines the variance of either type of curve. Must 
        have the same time units as `t`.
    tau_wt : float
        lifetime of decay of exponentially weighted Gaussian. Ignored if 
        `weighted` is set to False.
    tau : float
        lifetime of diffusion tail in the time units of `t`. Ignored if 
        `irf_type` is set to `'Gauss'`.
    b : float
        Relative contribution of diffusion tail to IRF. Ignored if 
        `irf_type` is set to `'Gauss'`.
    c : float
        Offset for IRF, recommended as zero.
    
    """
    if weighted == False:
        arg1 = Gauss(t, 1, t0, fwhm)
    else:
        arg1 = ExpGauss(t, t0, 1, fwhm, tau_wt, b=0)
    arg2 = np.exp(-(t-t0)/tau)
    for i, y in enumerate(arg2):
        if t[i] < t0:
            arg2[i] = 0 
    y = np.zeros(len(arg1))
    arg1 /= np.max(arg1)    
    for i in range(len(y)):
        y[i] = np.max([a*arg1[i], a*b*arg2[i]])
    c *=np.max(y)
    y[y<c] = c
    return y 
