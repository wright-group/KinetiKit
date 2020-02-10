"""
Functions for displaying things on the console
"""

import sys
import os
import csv

import numpy as np
import matplotlib.pyplot as plt

def printparamsexp(paramdict, errordict=None, decimals=4):
    """
    Prints the keys and values of a dictionary in an exponential form.
    
    Required Parameters
    -------------------
    paramdict : dictionary whose keys are the parameter names and values are
    the parameter values.
    paramdict : dictionary defining the errors of each parameter
    
    Optional Parameters
    -------------------
    decimals : integer
    How many decimal points to use to express the values. Default is 4. 
    
    """
    if errordict is None:
        for key, val in paramdict.items():
            print('{} = {:.{prec}e}'.format(key, val, prec=decimals))
    else:
        for key, val in paramdict.items():
            err = errordict[key]
            print('{} = {:.{prec}e} +/- {:.2e}'.format(key, val, err, prec=decimals))