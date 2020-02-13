#--- Structures on which different model types are built
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)

__all__ = ['RateModel', 'FunctionModel']

class RateModel:
    """
    Model structure based on rate equations
    """
    def __init__(self, **kwargs):
        params = self.default.copy()
        params.update(kwargs)
        
        self.name = self.class_name
        
        self.keys = params.keys()
        for key, value in params.items():
            setattr(self, key, value)

        self.popnum = len(self.populations)
    
    def rate(self, N, photons):
        raise NotImplementedError
        
    def PLsig(self, N):
        raise NotImplementedError
        
    def update(self, **kwargs):
        for key, val in kwargs.items():
            if key in self.keys:
                self.__setattr__(key,val)

    def params(self):
        return {key: self.__getattribute__(key) for key in self.keys}
    
    
class FunctionModel:
    """
    Model structure where signal is a function of time (e.g. multi-exponential)
    """
    def __init__(self, **kwargs):
        params = self.default.copy()
        params.update(kwargs)
        
        self.name = self.class_name
        
        self.keys = params.keys()
        for key, value in params.items():
            setattr(self, key, value)
            
        self.populations = None
        
    def PLsig(self, t):
        raise NotImplementedError
    
    def update(self, **kwargs):
        for key, val in kwargs.items():
            if key in self.keys:
                self.__setattr__(key,val)

    def params(self):
        return {key: self.__getattribute__(key) for key in self.keys}
        