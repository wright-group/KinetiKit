import numpy as np
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)

class MonoHTrap():
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
    
    def __init__(self, **kwargs):
        params = self.default.copy(); params.update(kwargs)
        self.name = self.class_name
        
        self.keys = params.keys()
        self.k_ann = params['k_ann'] # exciton annihilation
        self.k_dis = params['k_dis'] # exciton dissociation
        self.k_rec = params['k_rec'] # radiative free electron-hole recomb.fmono
        self.k_trx = params['k_trx'] # exciton trapping
        self.k_eth = params['k_eth'] # non-radiative el.-trapped hole rec.
        self.N_trp = params['N_trp'] # total number of traps
        self.cs = params['cs'] # absorption cross-section
        
        self.populations = ['x', 'e', 'h', 'th']
        self.popnum = len(self.populations)
   
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
    
    def update(self, **kwargs):
        for key, val in kwargs.items():
            if key in self.keys:
                self.__setattr__(key,val)

    def params(self):
        return {key: self.__getattribute__(key) for key in self.keys}
    
    
        return {key: self.__getattribute__(key) for key in self.keys}