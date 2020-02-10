import numpy as np
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)

# See instructions on how to program a model. Replace <text> with your desired parameters.

# A rough template:

# class yourModelName():
#     """
#     System description
#     """
    
#     class_name = 'yourModelName'
#     default = {'yourParameter1': 1,
#                'yourParameter2': 1,
#                'yourParameter3': 1,
#                'abs':1,
#                     }
    
#     def __init__(self, **kwargs):
#         params = self.default.copy(); params.update(kwargs)
#         self.name = self.class_name
        
#         self.keys = params.keys()
#         self.yourParameter1 = params['yourParameter1'] # parameter description
#         self.yourParameter2 = params['yourParameter2'] # parameter description
#         self.yourParameter3 = params['yourParameter3'] # parameter description
#         self.abs = params['abs'] # parameter description
        
#         self.populations = ['pop1', 'pop2', 'pop3'] # identify carrier types
#         self.popnum = len(self.populations)
   
#     def rate(self, N, photons):
#         pop1, pop2, pop3 = N 
        
#         pop1_rate = photons*self.abs - (self.yourParameter1)* pop1 #+ ...
#         pop2_rate = self.yourParameter1 * pop1 - self.yourParameter2 * pop2**2 #...
#         pop3_rate = #...
        
#         return np.array([pop1_rate, 
#                          pop2_rate, 
#                          pop3_rate])
        
#     def PLsig(self, N):
#         pop1 = N[self.populations.index('pop1')] 
        
#         return self.yourParameter1 * pop1 
           
#     def update(self, **kwargs):
#         for key, val in kwargs.items():
#             if key in self.keys:
#                 self.__setattr__(key,val)

#     def params(self):
#         return {key: self.__getattribute__(key) for key in self.keys}
