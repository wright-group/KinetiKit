import numpy as np
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)

class MonoEEA():
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
    
    def __init__(self, **kwargs):
        params = self.default.copy(); params.update(kwargs)
        self.name = self.class_name
        
        self.keys = params.keys()
        self.k_ann = params['k_ann'] # exciton annihilation
        self.k_dis = params['k_dis']
        self.k_eea = params['k_eea']
        self.k_rec = params['k_rec']
        self.cs = params['cs'] # absorption cross-section
        
        self.populations = ['x', 'e', 'h']
        self.popnum = len(self.populations)
   
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
           
    def update(self, **kwargs):
        for key, val in kwargs.items():
            if key in self.keys:
                self.__setattr__(key,val)

    def params(self):
        return {key: self.__getattribute__(key) for key in self.keys}
		
		
class MonoTransfer():
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
    
    def __init__(self, **kwargs):
        params = self.default.copy(); params.update(kwargs)
        self.name = self.class_name
        
        self.keys = params.keys()
        self.k_ann = params['k_ann']
        self.k_dis = params['k_dis']
        self.k_rec = params['k_rec']
        self.k_xtr = params['k_xtr']
        self.cs = params['cs']

        self.populations = ['x', 'e', 'h']
        self.popnum = len(self.populations)

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
        
    def update(self, **kwargs):
        for key, val in kwargs.items():
            if key in self.keys:
                self.__setattr__(key,val)

    def params(self):
        return {key: self.__getattribute__(key) for key in self.keys}
    