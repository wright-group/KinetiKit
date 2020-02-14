"""
You can use a custom-built model on any of the example files, by calling a 
system from a different file in the same working directory. You will need to 
explicitly import the models from their file:
"""

from custom_models import Triexp, MonoEEA
"""
If the above line does not work, make sure that the current working directory
is the directory where the custom_models.py file is located.
After importing the systems, you can call them by their name inside the file:
"""

system1 = Triexp()
print(system1.name, system1.keys)
system2 = MonoEEA()
print(system2.name, system2.keys)

