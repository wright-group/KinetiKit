import numpy as np
import warnings
from KinetiKit.sim.systems import RateModel
warnings.filterwarnings("ignore", category=RuntimeWarning)

class MonoHetero1(RateModel):
        """
        Model used for simultaneous fitting of a system while its pure
        and while it's in a heterostructure. Contributor: Jing Li 2020
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