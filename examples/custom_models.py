"""
Example file a user of the program would write to define new models.
To use such a file, you should place it in the same directory as the 
main file that you edit to import and simulate data.

"""

import numpy as np
import warnings
from KinetiKit.sim.systems import RateModel, FunctionModel

warnings.filterwarnings("ignore", category=RuntimeWarning) #use this to avoid annoying messages

class MonoEEA(RateModel):
    """
    Single (monolithic) system, excitons, free carriers, exciton-exciton 
    annihilation possible.
    """
    
    class_name = 'MonoEEA'
    default = {'k_ann': 1,
               'k_dis': 1,
               'k_eea': 1,
               'k_rec':1,
               'cs': 1,
                    }
    
    populations = ['x', 'e', 'h']
    
    def rate(self, N, photons):
        nx, ne, nh = N # free excitons, trapped excitons
        
        nx_rate = photons*self.cs - (self.k_ann + self.k_dis) * nx - 0.5 * self.k_eea * nx**2 \
        + self.k_rec*ne*nh
        ne_rate = self.k_dis * nx - self.k_rec * ne * nh
        nh_rate = ne_rate
        
        return np.array([nx_rate, 
                         ne_rate, 
                         nh_rate])
        
    def PLsig(self, N):
        nx = N[self.populations.index('x')] 
        
        return self.k_ann * nx 

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
