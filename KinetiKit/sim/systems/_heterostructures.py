import numpy as np
import warnings
from KinetiKit.sim.systems import RateModel
warnings.filterwarnings("ignore", category=RuntimeWarning)
    
class Hetero(RateModel):
    """
    Heterostructure of two adjacent systems, excitons, free carriers, and 
    transfer thereof. No trapping.
    """
    
    class_name = 'Hetero'
    default = {
    'k1_ann': 1,
    'k1_dis': 1,
    'k1_rec': 1,
    'cs1': 1,
    
    'k2_ann': 1,
    'k2_dis': 1,
    'k2_rec': 1,
    'cs2': 1,
    
    'k_xtr': 1,
    'k_etr': 1,
    'k_htr': 1,
    }
    
    populations = ['1x', '1e', '1h', '2x', '2e', '2h']
   
    def rate(self, N, photons):
        n1x, n1e, n1h, n2x, n2e, n2h = N
        n1x_rate = photons*self.cs1 \
            - (self.k1_ann + self.k1_dis + self.k_xtr) * n1x
        n1e_rate = self.k1_dis*n1x - self.k1_rec*n1e*n1h - self.k_etr*n1e
        n1h_rate = self.k1_dis*n1x - self.k1_rec*n1e*n1h - self.k_htr*n1h

        n2x_rate = photons*self.cs2 - (self.k2_ann + self.k2_dis) * n2x \
            + self.k_xtr*n1x
        n2e_rate = self.k2_dis*n2x - self.k2_rec*n2e*n2h + self.k_etr*n1e
        n2h_rate = self.k2_dis*n2x - self.k2_rec*n2e*n2h + self.k_htr*n1h

        return np.array([n1x_rate,
                         n1e_rate,
                         n1h_rate,
                         n2x_rate,
                         n2e_rate,
                         n2h_rate])

    def PLsig(self, N):
        # This produces two PL outputs - one for each part of the heterostructure
        out = np.zeros((2, N.shape[-1]))
        n1x = N[self.populations.index('1x')]
        n1e = N[self.populations.index('1e')]
        n1h = N[self.populations.index('1h')]
        n2x = N[self.populations.index('2x')]
        n2e = N[self.populations.index('2e')]
        n2h = N[self.populations.index('2h')]
        out[0] = self.k1_ann * n1x + self.k1_rec * n1e * n1h
        out[1] = self.k2_ann * n2x + self.k2_rec * n2e * n2h
        return out
    