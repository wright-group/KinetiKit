import numpy as np
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
from KinetiKit.sim.systems import RateModel

class Mono(RateModel):
    """
    Single (monolithic) system, excitons, free electrons & holes, no trapping.
    Electrons and holes recombine to produce light.
    """
    
    class_name = 'Mono'
    default = {'k_ann': 1,
               'k_dis': 1,
               'k_rec': 1,
               'cs': 1,
                    }
    populations = ['x', 'e', 'h']
    
    def rate(self, N, photons):
        nx, ne, nh = N
        nx_rate = photons*self.cs - (self.k_ann + self.k_dis) * nx
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


class MonoRecX(RateModel):
    """
    Same as Mono, but electrons and holes can recombine into excitons, instead 
    of just producing light.
    """
    
    class_name = 'MonoRecX'
    default = {
        'k_ann': 1,
        'k_dis': 1,
        'k_rec': 1,
        'cs': 1,
        } 
    
    populations = ['x', 'e', 'h']
    
    def rate(self, N, photons):
        nx, ne, nh = N
        nx_rate = photons*self.cs - (self.k_ann + self.k_dis) * nx + self.k_rec*ne*nh
        ne_rate = self.k_dis*nx - self.k_rec*ne*nh
        nh_rate = ne_rate

        return np.array([nx_rate,
                         ne_rate,
                         nh_rate])

    def PLsig(self, N):
        nx = N[self.populations.index('x')] 

        return self.k_ann * nx
    
class MonoFracFree(RateModel):
    """
    Same as MonoRecX, but now excitation also produces a fraction of free
	carriers in addition to excitons. Electrons and holes recombine into excitons.
    """
    
    class_name = 'MonoFracFree'
    default = {'k_ann': 1,
               'k_dis': 1,
               'k_rec': 1,
               'r_free': 1,
               'cs': 1,
                    }
    populations = ['x', 'e', 'h']
    
    def rate(self, N, photons):
        nx, ne, nh = N
        nx_rate = photons*self.cs*(1-self.r_free) - (self.k_ann + self.k_dis)*nx + self.k_rec*ne*nh
        ne_rate = photons*self.cs*self.r_free + self.k_dis*nx - self.k_rec*ne*nh
        nh_rate = ne_rate

        return np.array([nx_rate,
                         ne_rate,
                         nh_rate])

    def PLsig(self, N):
        nx = N[self.populations.index('x')] 

        return self.k_ann * nx



    
    
