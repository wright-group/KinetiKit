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
	
class MonoHetero1(RateModel):
        """
        Model used for simultaneous fitting of a system while its pure
        and while it's in a heterostructure. Contributor: Jing Li 2020
        Two PL outputs are considered: pure compound and same compound
        inside heterostructure.
        """
    
        class_name = 'MonoHetero1'
        default = {
            'k_ann': 1,
            'k_trp': 1,
            'k_dtr': 1,
            'k_xtr': 1,
            'cs': 1,
        }
        
        populations = ['xp', 'tp', 'xh', 'th']
            
        def rate(self, N, photons):
            nxp, ntp, nxh, nth = N  #carrier species
            nxp_rate = photons*self.cs - (self.k_ann + self.k_trp) * nxp+self.k_dtr*ntp
            ntp_rate = -self.k_dtr*ntp + self.k_trp*nxp
            nxh_rate = photons*self.cs - (self.k_ann + self.k_trp) * nxh+self.k_dtr*nth-self.k_xtr*nxh
            nth_rate = -self.k_dtr*nth + self.k_trp*nxh
    
            return np.array([nxp_rate,
                             ntp_rate,
                             nxh_rate,
                             nth_rate])
    
        def PLsig(self, N):
            out = np.zeros((2, N.shape[-1]))
            nxp = N[self.populations.index('xp')] 
            nxh = N[self.populations.index('xh')]
            out[0] = self.k_ann * nxp
            out[1] = self.k_ann * nxh
    
            return out
	


