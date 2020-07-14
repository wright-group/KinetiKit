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
    Data object created by importing one time-resolved .asc signal data file from
    a Becker & Hickl SPC-130 single photon counting module (SPCM software)
    
    Can be used to represent a generic TRPL text file import, but the 
    `metadata` keyword cannot be used meaningfully in this case.
    
    Parameters
    ----------
    filename : string
        the file path of interest
    
    Optional Parameters
    -------------------
    key : string
        an identifier of choice for the object, which can be later used in 
        labels etc. Default is `None`.
    time_unit : string from ``KinetiKit.units.units`` dictionary keys
        The unit in which time is expressed in the data file. Default is`'ns'`.
    coll : float
        Collection time of the experiment in seconds. If `metadata` is `True`, 
        and the specific Becker & Hickl SPCM format is used,`coll` is 
        over-ridden by the collection time embedded in the file. Default is 1.
    weigh_by_coll : boolean
        If `True`, the signal values (y) of the object are automatically 
        divided by the collection time to yield counts-per-second units. 
        Default is `True`.
    metadata : boolean
        If `True` the collection time is obtained from inside the .asc file.
        For this data object,`metadata=True` only works for Becker & Hickl SPCM
        data files. Default is  `False`.
    pulse_power : float
        pulse power information (in any desired units) for the data. Default is
        `None`.
    cw_power : float
        continuous-wave power information (in any desired units) for the data.
        Default is `None`.
    wavelength : float
        wavelength information (in any desired units) for the data. Default is
        `None`.
    delimiter : string
        Delimiter in the text file. Default is `','`.
    unpack : boolean
        See `unpack` argument of ``numpy.genfromtxt``.
    skip_h : integer
        See `skip_header` argument of ``numpy.genfromtxt``.
    skip_f : integer
        See `skip_footer` argument of ``numpy.genfromtxt``.
        
    """
    
    def __init__(self, filename, key=None, time_unit='ns', coll=1,
                 metadata=False, weigh_by_coll=True, 
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
        """
        Wrapper function for the ``scipy.signal.savgol_filter`` smoothing 
        function with limited relevant arguments
        
        """
                                
        self.y = signal.savgol_filter(self.y, window_length = window_length, 
                                    polyorder = polyorder, mode = mode)
    
    def interp(self, t):
        """
        Interpolates the y-values of the imported object along a time axis.
        Interpolation is most accurate if the first and last element of the 
        time array match the first and last element of the object's x-data.
        
        Parameters
        ----------
        t : 1-D array or dictionary of the type ``sim.time.linear()``
            Defines the time axis along which to interpolate.

        """
        if isinstance(t, numpy.ndarray):
            time_array = t
            self.y = np.interp(time_array, self.x, self.y)
            self.x = time_array
            
        if isinstance(t, dict):
            time_array = t['array'][::t['subsample']]
            self.y = np.interp(time_array, self.x, self.y)
            self.x = time_array
        
    
    def remove_zeros(self):
        """
        Replaces any negative or zero elements of the object's y-data with the 
        minimum nonzero positive value of the y-data array.
        """
        
        y=self.y
        ymin = min(y[np.where(y>0)])
        ynew = y.copy()
        ynew[np.where(ynew<=0)] = ymin
        self.y = ynew
    
    def extract_average(self, t_start=12*ns, t_end=12.5*ns):
        """
        Returns the average y-value of the object defined between two x-values.
        The x-values should be given in seconds and `t_start` must be less than
        `t_end`.
        
        """
        
        i_start, t_start = kin_kit.find_nearest(self.x, t_start)
        i_end, t_end = kin_kit.find_nearest(self.x, t_end)
        return np.average(self.y[i_start:i_end])
        
    def extract_peak(self, avgnum=5):
        """
        Returns the peak y-value of the object, averaging between a defined 
        number of points (`avgnum`) to control for noise.
        
        """
        peak_values = []
        y = list(self.y)
        for i in range(avgnum):
            peak_values.append(max(y))
            y.remove(max(y))
        return np.average(peak_values)
    
    def copy(self):
        """
        Returns an independent copy of the data object and all its attributes
        at the time this function is called.
        """
        return copy.deepcopy(self)
        
        
