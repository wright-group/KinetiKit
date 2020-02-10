"""
Methods for reading csv files created by this program.
"""

import sys
import os
import csv

import numpy as np

def csv_to_dict(filename, sourcefolder='', headertype='row'):
    """
    Returns a dictionary with keys and values obtained from a CSV file. Useful 
    for importing fit parameters and X-Y traces from previous fits.
    
    Parameters
    ----------
    filename : string
        Filename of CSV file containing the data and headings to be imported.
    sourcefolder : path (or string)
        Location of CSV file.
    headertype : 'row' or 'column'
        Specifies how the headings are arranged in the CSV file. In a
        '*_params.csv' file, each heading takes up a `row`. In a '*_traces.csv'
        file, each heading takes up a column, with arrays following underneath.
    """
    
    p = os.path.join(sourcefolder, filename)
    
    if headertype == 'row':
        # assumes that the headers are arranged in each row, with a single 
        # value following to the right of each header (e.g. a _params file)
        headers = np.genfromtxt(p, delimiter = ',', unpack=True, usecols=[0],
                                dtype=str) 
        values = np.genfromtxt(p, delimiter = ',', unpack=True, usecols=[1],
                               dtype=float)
    if headertype == 'column':
        # assumes that the headers are arranged in each column, with a single
        # value or an array of values stacked vertically below each header
        with open(p, 'r') as f:
            reader = csv.reader(f, delimiter=',')
            # get header from first row
            headers = next(reader)
            # get all the rows as a list
            data = list(reader)
            # transform data into numpy array
            values = np.transpose(np.array(data).astype(float))
    
    outdict = {}
    
    for i, header in enumerate(headers):
        outdict[header] = values[i]
    
    return outdict
