#--- Units
# The letter 'u' here is meant to replicate the Greek letter Î¼ (mu)
# If you absolutely hate this notation, alternatives using 'micro' instead of u
# are provided

fs = 1e-15
ps = 1e-12
ns = 1e-9
us = 1e-6
microsec = 1e-6
ms = 1e-3
MHz = 1e6
kHz = 1e3
nW = 1e-9
uW = 1e-6
microWatt = 1e-6
nm = 1e-9
um = 1e-6
micron = 1e-6
mm = 1e-3
cm = 1e-2
eV = 1.602e-19
meV = 1.602e-22

#--- Unit Dictionary (useful as an input on functions)
units = {'fs' : fs,
         'ps' : ps,
         'ns' : ns,
         'us' : us,
         'microsec' : microsec,
         'μs' : microsec,
         'ms' : ms,
         'MHz' : MHz, 
         'kHz' : kHz,
         'nm' : nm,
         'μm' : um,
         'nW' : nW,
         'uW' : uW,
         'μW' : uW,
         'microWatt' : microWatt,
         'micron' : micron,
         'mm' : mm,
         'cm' : cm,
         'eV' : eV,
         'meV': meV,
         'none' : 1,
        }

#--- constants
c = 2.9979e8 # speed of light in m/s
hPlanck = 6.626e-34 # Planck's constant in J s
k_B = 1.3806e-23 # Boltzmann's constant in J/K
