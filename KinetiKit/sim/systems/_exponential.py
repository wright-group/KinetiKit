import numpy as np
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)

class Biexp():
    """
    Heterostructure of two adjacent systems, excitons, free carriers, and 
    transfer thereof. No trapping. To be used with sim.lib.simulate_func 
    instead of sim.lib.refined_simulation.
    """
    
    class_name = 'Biexp'
    default = {
    'A1': 1,
    'tau1': 1,
    'tau2': 1,
    'offset':0, 
    }

    def __init__(self, **kwargs):
        params = self.default.copy(); params.update(kwargs)
        self.name = self.class_name
        self.keys = params.keys()
        
        self.A1 = params['A1']
        self.tau1 = params['tau1']
        self.tau2 = params['tau2']
        self.offset = params['offset']
        
        self.populations = None # need to define this for non-kinetic models as
                                # it determines how `simulate_and_compare` 
                                # treats the system.
        
    def PLsig(self, t):
        return self.A1*np.exp(-t/self.tau1) + (1-self.A1)*np.exp(-t/self.tau2) + self.offset
    
    def update(self, **kwargs):
        for key, val in kwargs.items():
            if key in self.keys:
                self.__setattr__(key,val)
    
    def params(self):
        return {key: self.__getattribute__(key) for key in self.keys}
