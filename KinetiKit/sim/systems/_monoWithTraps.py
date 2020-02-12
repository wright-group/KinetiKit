import numpy as np
import warnings
from KinetiKit.sim.systems import RateModel
warnings.filterwarnings("ignore", category=RuntimeWarning)

class MonoHTrap(RateModel):
    """
    Single (monolithic) system, excitons, excitons trapped as holes.
    """
    
    class_name = 'MonoHTrap'
    default = {
        'k_ann': 1,
        'k_dis': 1,
        'k_rec': 1,
        'k_trx': 1,
        'k_eth': 1,
        'N_trp': 1,
        'cs': 1,
        }  
    
    populations = ['x', 'e', 'h', 'th']
   
    def rate(self, N, photons):
        nx, ne, nh, nth = N # free excitons, trapped excitons
        
        nx_rate = photons*self.cs - self.k_ann * nx - self.k_dis * nx \
        - self.k_trx * (self.N_trp - nth) * nx + self.k_rec * ne * nh
        
        ne_rate = self.k_dis * nx - self.k_rec * ne * nh \
        + self.k_trx * (self.N_trp - nth) * nx - self.k_eth * nth * ne 
        
        nh_rate = self.k_dis * nx - self.k_rec * ne * nh
        
        nth_rate = self.k_trx * (self.N_trp - nth) * nx - self.k_eth * nth * ne
        
        return np.array([nx_rate,
                         ne_rate,
                         nh_rate,
                         nth_rate])
        
    def PLsig(self, N):
        nx = N[self.populations.index('x')] 
        
        return self.k_ann * nx 
