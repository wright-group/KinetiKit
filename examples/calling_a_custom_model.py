"""
You can use a custom-built model on any of the example files, by including the 
import process below, and, when defining the "system" in your code, replacing 
sim.systems.<modelname> with myModels.<modelname>.

Example folder structure for this to work:
    
SomeFolder
|-- calling_a_custom_model.py
|-- UserDefinedModels
|   |-- custom_models.py
|   |-- __init__.py
|   |  | from custom_models import *

"""

import UserDefinedModels as myModels

system = myModels.Triexp()
system.update(**{'A1':0.9, 'tau1': 1e10, 'tau2':1e8, 'offset':0.2})
print(system.params())