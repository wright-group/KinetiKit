# -*- coding: utf-8 -*-
"""
Functions that output an overlay of simulation and data given a system, 
a set of y data, and parameters. 

@author: Natalia Spitha
"""

import numpy as np
from ipywidgets import interact, fixed
import ipywidgets as widgets

from KinetiKit import kit as kin_kit
from KinetiKit.units import units, ns

def active_widgets(system, dictionary, total_args=20):
    """
    This function:
        - accepts a dictionary of widgets each corresponding to a 
        system parameters
        - ensures the dictionary has the parameters in the correct order, 
        i.e. the order defined in sim.systems.ChosenSystem.
        - returns a list of widgets, with size total_args. This list contains
        the widgets for the system parameters to be displayed interactively, 
        and contains inactive (fixed(0)) widgets for the rest of the slots 
        in the list.
        
    """
    
    if system.params().keys() != dictionary.keys():
        print('Ensure that \'parameters\' contains ALL parameters of the system.')
        return
    
    else:
        ordered_dict = system.default.copy()
        ordered_dict.update(dictionary)
        dictionary = ordered_dict
        
        all_widgets = list(np.zeros(total_args))
        param_list = list(dictionary.keys())
        for i, widg in enumerate(all_widgets):
            if i < len(param_list):
                all_widgets[i] = dictionary[param_list[i]]
            else:
                all_widgets[i] = fixed(0)
        # all_widgets[np.index(all_widgets==0] = fixed(0)
        return all_widgets

def powerwidget(slidepower=False, pmin=0.040, pmax=1.5, pset = 0.5, step=0.01, power=0, unit='μW'):
    if not slidepower:
        return fixed(power)
    else:
        return widgets.FloatSlider(value=pset, min=pmin, max=pmax, step=step, description = 'power (%s)'%unit,
                                   disabled=not slidepower, continuous_update=False, orientation='Horizontal',
                                   readout=slidepower, readoutformat='.3f')

def makeWidget(key, val, slider='linear', min_val=None, max_val=None, step=0.01, disabled=False,
               continuous_update=False, orientation='horizontal', readout=None,
               readout_format=None):
    
    #slide should be 'linear', 'log', or 'fixed'
    if readout is None:
        readout = not disabled
    
    if slider == 'linear':
        sliderfunc = widgets.FloatSlider
        if min_val is None:
            min_val = 0
        if max_val is None:
            max_val = 2*val
        if readout_format is None:
            precision = kin_kit.precision(step)
            readout_format = '.{pr}f'.format(pr=precision)
        if disabled:
            precision = kin_kit.precision(val)
            readout_format = '.{pr}f'.format(pr=precision)
            
            
    elif slider == 'log':
        sliderfunc = widgets.FloatLogSlider
        # FloatLogSlider accepts the actual value as `value`, but accepts 
        # log limits instead of actual power limits. So, the min_val and
        # max_val provided are converted into their log values before the
        # slider function is called.
        
        if min_val is None:
            min_val = np.log10(val)-3
        else: 
            min_val = np.log10(min_val) 
            
        if max_val is None:
            max_val = np.log10(val)+3
        else:
            max_val = np.log10(max_val)
        
        if readout_format is None:
            readout_format = '0.3e'
            
    elif slider == 'fixed':
        return fixed(val)
    
    return sliderfunc(value=val,
                          min=min_val,
                          max=max_val,
                          step=step,
                          description=key, 
                          disabled=disabled,
                          continuous_update=continuous_update,
                          orientation=orientation,
                          readout=readout,
                          readout_format=readout_format)
    
        