"""
Methods used to compile fits, traces, parameters, etc. for saving
"""

import sys
import os
import csv

import numpy as np
import matplotlib.pyplot as plt

def saveparam_dict(system, opt_DE, opt_LS, bounds, fit_args, doFit):
    """
    Creates a dictionary that contains the fitting results and parameters for
    a given fitting routine. Useful for saving those parameters using 
    ``save_fit``. Currently designed for to accept `opt` as output by the 
    scipy.differential_evolution algorithm.
    
    Parameters
    ----------
    system: a system instance
    opt_DE : dictionary or None
        Optimization Result from Differntial Evolution algorithm
    opt_LS : tuple or None
        Optimization Result from Least Squares algorithm
    bounds: dictionary
        Dictionary whose keys are the names parameters with respect to which 
        the optimization was performed. Each value of the dictionary should be
        a tuple of size 2 representing the lower and upper boundaries passed
        to the ``differential_evolution`` method.
    fit_args: tuple
        Tuple containing the arguments passed into the cost function. Designed
        for use with the output of ``fit.lib.sac_args``.
    doFit: Boolean
        Indicates whether the current parameters of the system were obtained
        via an optimization or not.     
    
    """
    
    saveparam_dict = system.params()
    saveparam_dict['Model type'] = system.name
    
    if doFit:
        if opt_LS is not None:
            saveparam_dict['Sum of Sq. of Residuals'] = np.sum(opt_LS[2]['fvec']**2)
        saveparam_dict['No. of DE Evaluations'] = opt_DE.nfev
        saveparam_dict['No. of DE Iterations'] = opt_DE.nit
        if opt_LS is not None:
            saveparam_dict['No. of LS Evaluations'] = opt_LS[2]['nfev']
        for i, bkey in enumerate(list(bounds.keys())):
            saveparam_dict[bkey+'_bounds'] = bounds[bkey]
        for i, bkey in enumerate(list(bounds.keys())):
            if opt_LS is not None:
                saveparam_dict[bkey+'_error'] = opt_LS[1][i]
        saveparam_dict['comparison'] = fit_args[-7]
        saveparam_dict['absolute'] = fit_args[-6]
        saveparam_dict['time_limits'] = fit_args[-5]
        saveparam_dict['norm_comp'] = fit_args[-4]
        saveparam_dict['roll_criterion'] = fit_args[-3]
        
    return saveparam_dict

def dict_to_csv(filename, dictionary, paramdict=True):
    """
    Writes a .csv file listing the entries from a dictionary.
    
    Parameters
    ----------
    filename : ``string``
        Containing the desired filename WITHOUT the .csv extension
    dictionary : dictionary
        Will be saved item-by-item in a csv file
    paramdict : boolean, optional
        A statement asking whether the dictionary to be converted contains
        fitting parameters, in which case the string '_params' is appended to
        `filename`. Default is true.
    
    """
    if paramdict:
        filename += '_params'
    
    with open(filename + '.csv', 'w') as f:
        for key in dictionary.keys():
            if type(dictionary[key]) == tuple:
                f.write("%s,%s, %s \n"%(key, dictionary[key][0], dictionary[key][1]))
            else:
                f.write("%s,%s \n"%(key, dictionary[key]))
        print('Saved %s.csv.'%filename)
   
def traces_to_csv(filename, dictionary):
    """
    Writes a .csv file containing x and y data from a dictionary.
    
    Parameters
    ----------
    filename : ``string``
        Containing the desired filename WITHOUT the .csv extension
    dictionary : dictionary
        The keys of the dictionary must contain the title of each trace
        (or the time), and the ditionary values must be arrays of the same
        size.
    """
    
    headings = []
    for key_string in list(dictionary.keys()):
        headings.append(key_string)
    
    # Use first element of dict to determine array length
    length = len(dictionary[headings[0]])
    filename += '_traces'
    with open(filename + '.csv', 'w') as f:
        f.write(','.join(headings))
        f.write('\n')
        for i in range(length):
            values = []
            for key in dictionary.keys():
                values.append(dictionary[key][i].astype(str))
            f.write(','.join(values))
            f.write('\n')
        print('Saved %s.csv.'%filename)

def save_fit(filename, traces={}, params={}, destinationfolder=''):
    """
    A shortcut for performing ``traces_to_csv`` and ``dict_to_csv``
    to save the fit parameters and data/sim traces of successful fits, and
    saving any necessary figures. 
    The function will add the number 1 to the end of the provided filename,
    and increase that number if the filename already exists
    """
    
    od = os.path.join(destinationfolder, filename) + '_fit'
    
    for i in np.arange(0,500,1):
        d = od + str(i+1)
        print(d)
        
        condition = os.path.exists(d+'.svg') or os.path.exists(d+'.png') or os.path.exists(d+'_params.csv') or os.path.exists(d+'_traces.csv') 
        if condition:
            pass
        
        else:
            dict_to_csv(d,params)
            traces_to_csv(d,traces) 
            plt.savefig(d+'.svg', format='svg')
            print('Saved ' + d + '.svg')
            plt.savefig(d+'.png', format='png')
            print('Saved ' + d + '.png')
            break

def save_species(filename, destinationfolder='', norm=False):
    """
    Currently, this function saves only the last-generated figure.  The 
    function will add the number 1 to the end of the provided filename,
    and increase that number if the filename already exists
    """
    
    od = os.path.join(destinationfolder, filename) + '_fit'
    
    for i in np.arange(0,500,1):
        d = od + str(i+1)
        if norm:
                d = d + '_norm'
        condition = os.path.exists(d+'_species.svg') or os.path.exists(d+'_species.png')
        if condition:
            pass
			
        else:
            plt.savefig(d+'_species.svg', format='svg')
            print('Saved ' + d + '_species.svg')
            plt.savefig(d+'_species.png', format='png')
            print('Saved ' + d + '_species.png')
            break

    
