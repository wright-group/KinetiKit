"""
Test file for fitting a data object with a simulation of the Mono type. 
"""

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from KinetiKit import sim, fit, data, artists
from KinetiKit import kit as kin_kit
from KinetiKit.units import nW, uW, MHz, nm, ps, ns

#--- Creating Time Object
to = sim.time.linear(N=1000, period=12.5*ns)
dtime = to['array'][::to['subsample']]

#--- Create System Instance (select model from sim.systems)
system = sim.systems.Mono()

#--- Parameters of simulation
params = {
    'k_ann': 1.25e8,
    'k_dis': 2.01e9,
    'k_rec': 1.3e3,
    'cs': 0.5,
    }

power = 1000 * nW
system.update(**params)

#--- Create Excitation object
light = sim.lib.Excitation(pulse={'power':power,
                                  'reprate': 80*MHz,
                                  'wavelength': 400*nm})

transient, converged = sim.lib.refined_simulation(system, to, light) # system re-simulated with fit parameters
species_set = transient # set of population's time evolution
pl = system.PLsig(transient) # PLsig function applied to population arrays
sims = sim.lib.convolve_irf(pl, dtime, fwhm=40*ps) # PL signal convolved with IRF

# Aligns sim either by max. (align_by_max) or steep (align_by_steep)
align_to = 0.5*ns #value at which data and simulation are aligned for saving/diplay

aligned_sims = kin_kit.align_by_max(sims, dtime, value = align_to, avgnum=5)

#--- Plot 
artists.plot3scales.plot(dtime, aligned_sims, system, t_dict=to, ivtype='none', 
                     annotate=True, fig=None, ResetColorCyc=True, 
                     linewidth=6, opaq=0.4,
                     mlabel='mono sim')
   
def plot_species(save=False):
    artists.plot3scales.show_species(dtime, [species_set], system, False, 
                                     save=save)
