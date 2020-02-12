import numpy as np
import warnings

from KinetiKit.sim.systems import FunctionModel

warnings.filterwarnings("ignore", category=RuntimeWarning)

class Biexp(FunctionModel):
    """
    Biexponential fitting. To be used with sim.lib.simulate_func 
    instead of sim.lib.refined_simulation. 
    """
    
    class_name = 'Biexp'
    default = {
    'A1': 1,    # Amplitude of first lifetime
    'tau1': 1,  # First lifetime
    'tau2': 1,  # Second lifetime
    'offset':0, # y-offset
    }
        
    def PLsig(self, t):
        return self.A1*np.exp(-t/self.tau1) + (1-self.A1)*np.exp(-t/self.tau2) \
            + self.offset

class Triexp(FunctionModel):
    """
    Triexponential fitting. To be used with sim.lib.simulate_func 
    instead of sim.lib.refined_simulation. 
    """
    
    class_name = 'Triexp'
    default = {
    'A1': 1,
    'tau1': 1,
    'A2' : 1,
    'tau2': 1,
    'tau3': 1,
    'offset': 1, 
    }
        
    def PLsig(self, t):
        return self.A1*np.exp(-t/self.tau1) + \
            self.A2*np.exp(-t/self.tau2) + \
                (1-self.A1 - self.A2) * np.exp(-t/self.tau3) + \
                    self.offset