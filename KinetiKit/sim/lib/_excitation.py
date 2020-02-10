import time

import numpy as np
import scipy
from scipy.signal import convolve

from KinetiKit import sim
from KinetiKit.units import units, MHz, fs, nm, uW, ps
import KinetiKit.kit as kin_kit


__all__ = ['Excitation']


class Excitation():
    """define excitation parameters for simulation"""
    
    pulse_default = {'power': 1*uW,
               'reprate': 80*MHz,
               'fwhm': 100*fs,
               'wavelength': 400*nm,
               'pulse_window': None,
               }
    
    cw_default = {'power': 0,
               'wavelength': 514.5*nm
               }
    
    def __init__(self, pulse={}, cw={}, numcycles=500):
        
        pulse_params = self.pulse_default.copy(); pulse_params.update(pulse)
        self.pulse = pulse_params
        self.pulse_power = self.pulse['power'] 
        self.pulse_energy = self.pulse['power'] / self.pulse['reprate']
        self.pulse_wavelength = self.pulse['wavelength']
        self.pulse_carriers = self.pulse_energy \
            * joules_to_photons(self.pulse_wavelength)
        self.pulse_fwhm = self.pulse['fwhm']
        self.pulse_window = self.pulse['pulse_window']
        if self.pulse_window is None:
            self.pulse_window = 20 * self.pulse_fwhm
        
        cw_params = self.cw_default.copy(); cw_params.update(cw)
        self.cw = cw_params
        self.cw_wavelength = self.cw['wavelength']
        self.cw_power = self.cw['power'] * joules_to_photons(self.cw_wavelength)
        
        self.numcycles = numcycles
        
    def gen_pulse(self, t, center=0, stepsize=None):
        pulseCarriers = self.pulse_carriers
        sigma_to_fwhm = 2 * np.sqrt(2 * np.log(2))
        pulsePeak = pulseCarriers * sigma_to_fwhm \
            / (np.sqrt(2 * np.pi) * self.pulse_fwhm)
        if len(t) > 50:
            return kin_kit.Gauss(t, pulsePeak, center, self.pulse_fwhm)
        else:
            return np.array([pulseCarriers/stepsize])
    
    def updated_with(self, pulse={}, cw={}, numcycles=None):
        """returns a new light object based on the current light object, but
        with some modified arguments.
        
        """
        new_pulse = self.pulse; new_pulse.update(pulse)
        new_cw = self.cw; new_cw.update(cw)
        if numcycles is None:
            numcycles = self.numcycles
        return Excitation(new_pulse, new_cw, numcycles)
        

def joules_to_photons(wavelength):
    """convert Joules to photons

    arguments
    wavelength : units of meters
    """
    h_Planck = 6.626e-34  # J s
    c = 2.9979e8  # m / s
    return wavelength / (h_Planck * c)
