import numpy as np
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)

class Mono():
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
    
    def __init__(self, **kwargs):
        params = self.default.copy(); params.update(kwargs)
        self.name = self.class_name
        
        self.keys = params.keys()
        self.k_ann = params['k_ann']
        self.k_dis = params['k_dis']
        self.k_rec = params['k_rec']
        self.cs = params['cs']

        self.populations = ['x', 'e', 'h']
        self.popnum = len(self.populations)
    
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
           
    def update(self, **kwargs):
        for key, val in kwargs.items():
            if key in self.keys:
                self.__setattr__(key,val)

    def params(self):
        return {key: self.__getattribute__(key) for key in self.keys}


class MonoRecX():
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
    
    def __init__(self, **kwargs):
        params = self.default.copy(); params.update(kwargs)
        self.name = self.class_name
        
        self.keys = params.keys()
        self.k_ann = params['k_ann'] # exciton annihilation
        self.k_dis = params['k_dis'] # exciton dissociation into e, h
        self.k_rec = params['k_rec'] # electron-hole recombination to light
        self.cs = params['cs'] # absorption cross-section
        
        self.populations = ['x', 'e', 'h']
        self.popnum = len(self.populations)
    
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
        ne = N[self.populations.index('e')] 
        nh = N[self.populations.index('h')]

        return self.k_ann * nx
        
    def update(self, **kwargs):
        for key, val in kwargs.items():
            if key in self.keys:
                self.__setattr__(key,val)

    def params(self):
        return {key: self.__getattribute__(key) for key in self.keys}
    
class MonoFracFree():
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
    
    def __init__(self, **kwargs):
        params = self.default.copy(); params.update(kwargs)
        self.name = self.class_name
        
        self.keys = params.keys()
        self.k_ann = params['k_ann']
        self.k_dis = params['k_dis']
        self.k_rec = params['k_rec']
        self.r_free = params['r_free']
        self.cs = params['cs']

        self.populations = ['x', 'e', 'h']
        self.popnum = len(self.populations)
    
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
        ne = N[self.populations.index('e')]
        nh = N[self.populations.index('h')]

        return self.k_ann * nx# + self.k_rec * ne * nh
           
    def update(self, **kwargs):
        for key, val in kwargs.items():
            if key in self.keys:
                self.__setattr__(key,val)

    def params(self):
        return {key: self.__getattribute__(key) for key in self.keys}


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
    
    
    
