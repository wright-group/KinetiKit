"""
Functions and class definitions for handling experimental data.
Author: Natalia Spitha, 6/3/2019
"""

import numpy as np
from scipy import signal

from KinetiKit import sim
from KinetiKit.units import units, ns
from KinetiKit import kit as kin_kit
import warnings
import numpy
import copy

__all__ = ['data_from_SPCM']


class data_from_SPCM(object):
    """
    Data object created by importing one time-resolved signal data file from
    a Becker & Hickl SPC-130 single photon counting module. 
    
    Can be used to represent a generic TRPL text file import, but the 
    `metadata` keyword cannot be used meaningfully in this case.
    """
    
    def __init__(self, filename, key=None, time_unit='ns', coll=1,
                 weigh_by_coll=True, metadata=False,
                 pulse_power=None, cw_power=None, wavelength=None,
                 delimiter=',', unpack=True, skip_h=10, skip_f=1):
        
        self.filename = filename
        self.time_unit = time_unit
        self.key = key # identifier
        self.weigh_by_coll = weigh_by_coll
        self.skip_h = skip_h; self.skip_f = skip_f
        
        # initializing some attributes
        self.dark_subtracted=False
        self.weighed_by_coll = False
        
        # If aquisition information was included in file, this part extracts
        # the aquisition time:
        if metadata:
            with open(self.filename) as file:
                line_col = file.readlines()[45]
                self.coll = float(line_col.split(',')[2].split(']')[0])
    
        else:
            self.coll = coll
                
        self.x, self.y = np.genfromtxt(filename, delimiter=',', unpack=unpack, skip_header=skip_h, skip_footer=skip_f) 
        
        if weigh_by_coll:
            self.y /= self.coll
            self.weighed_by_coll =True
                
        self.x -= self.x[0]
        self.x *= units[time_unit]
        
        self.pulse_power = pulse_power
        self.cw_power = cw_power
        
        self.wavelength = wavelength
        
    
    def dark_subtract(self, dark_counts, method='average'):
        """
        Accepts a float, array, or data object representing dark counts, and
        performs a dark count subtraction on the current data object.
        
        This function does not check for intensity or time units; ensure that
        the count units (recommended: Counts Per Second) are consistent between
        the data object and the dark count instance.
        
        Parameters
        -------------------
        dark_counts : float, array, or KinetiKit.data.lib.data object
            If dark_counts is an array, it should have the same length as 
            the y array of the data object.
        method : 'average' or 'elementwise', optional
        
        """
        
        if not self.dark_subtracted:
            
            if isinstance(dark_counts, numpy.float) or isinstance(dark_counts, numpy.int):
                self.y -= dark_counts
            elif isinstance(dark_counts, numpy.ndarray):
                if method == 'elementwise': 
                    self.y -= dark_counts
                elif method == 'average':
                    self.y -= np.average(dark_counts)
            elif isinstance(dark_counts, data_from_SPCM):
                if self.weighed_by_coll==False or dark_counts.weighed_by_coll==False:
                    warnings.warn('Dark subtraction was performed. However, we recommend that you convert Counts to Counts Per Second by weighing data objects by the collection time.', UserWarning)
                if method == 'elementwise':
                    self.y -= dark_counts.y
                elif method == 'average':
                    self.y -= np.average(dark_counts.y)
            else:
                return
            self.dark_subtracted=True
            
        else: print('Dark subtraction has already been performed')
        
        
    def smooth(self, window_length=17, polyorder=2, mode='nearest'):
        self.y = signal.savgol_filter(self.y, window_length = window_length, 
                                    polyorder = polyorder, mode = mode)
    
    def interp(self, t):
        if isinstance(t, numpy.ndarray):
            time_array = t
            self.y = np.interp(time_array, self.x, self.y)
            self.x = time_array
            
        if isinstance(t, dict):
            time_array = t['array'][::t['subsample']]
            self.y = np.interp(time_array, self.x, self.y)
            self.x = time_array
        
    
    def remove_zeros(self):
        y=self.y
        ymin = min(y[np.where(y>0)])
        ynew = y.copy()
        ynew[np.where(ynew<=0)] = ymin
        self.y = ynew
    
    def extract_average(self, t_start=12*ns, t_end=12.5*ns):

        y = self.y
        i_start, t_start = kin_kit.find_nearest(self.x, t_start)
        i_end, t_end = kin_kit.find_nearest(self.x, t_end)
        return np.average(self.y[i_start:i_end])
        
    def extract_peak(self, avgnum=5):
        peak_values = []
        y = list(self.y)
        for i in range(avgnum):
            peak_values.append(max(y))
            y.remove(max(y))
        return np.average(peak_values)
    
    def copy(self):
        return copy.deepcopy(self)
        
        
