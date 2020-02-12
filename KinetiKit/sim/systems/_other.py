import numpy as np
import warnings
from KinetiKit.sim.systems import RateModel

warnings.filterwarnings("ignore", category=RuntimeWarning)

class MonoEEA(RateModel):
    """
    Single (monolithic) system, excitons, NO free carriers, exciton-exciton 
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
		
		
class MonoTransfer(RateModel):
    """
    Single (monolithic) system, excitons, free electrons & holes, no trapping, 
    transfer to or from the system. Unlike Hetero(), only one system's output
	is considered.
    """
    
    class_name = 'MonoTransfer'
    default = {
    'k_ann': 1,
    'k_dis': 1,
    'k_rec': 1,
    'k_xtr': 1,
    'cs': 1,
    }
    
    populations = ['x', 'e', 'h']

    def rate(self, N, photons):
        nx, ne, nh = N
        nx_rate = photons*self.cs - (self.k_ann + self.k_dis) * nx - self.k_xtr * nx
        ne_rate = self.k_dis*nx - self.k_rec*ne*nh
        nh_rate = ne_rate

        return np.array([nx_rate,
                         ne_rate,
                         nh_rate])

    def PLsig(self, N):
        nx = N[self.populations.index('x')] 
        ne = N[self.populations.index('e')]
        nh = N[self.populations.index('h')]

        return self.k_ann * nx + self.k_rec * ne * nh
