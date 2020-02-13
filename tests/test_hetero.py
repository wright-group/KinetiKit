"""
Test file for fitting multiple data traces with a simulation of Hetero type, 
e.g. in a heterostructure with different emitting colors. 
Requires numpy, scipy, matplotlib, and time packages along with TRPL package

"""
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from KinetiKit import sim, fit, data, artists
from KinetiKit import kit as kin_kit
from KinetiKit.units import nW, uW, MHz, nm, ps, ns

#--- Creating Time Object
to = sim.time.linear(N=1000, period = 12.5*ns)
dtime = to['array'][::to['subsample']]

#--- Create System Instance (select model from sim.systems)
system = sim.systems.Hetero()

#--- Parameters of simulation 

initparams = {
        'k1_ann': 2.298e9,
        'k1_dis': 2.132e9,
        'k1_rec': 6.181e5,
        'cs1': 0.2,
        'k2_ann': 7.242e8,
        'k2_dis': 4.816e9,
        'k2_rec': 4.574e4,
        'cs2': 0.2,
        'k_xtr': 2.407e9,
        'k_etr': 0,
        'k_htr': 0,
    }  

power = 1000*nW
system.update(**initparams)

#--- Create Excitation object and initial system
light = sim.lib.Excitation(pulse={'power':power,
                                  'reprate': 80*MHz,
                                  'wavelength': 400*nm})
# Simulation
transient, converged = sim.lib.refined_simulation(system, to, light, N_coarse=500)
species_set = transient
pl = system.PLsig(transient)
sims = sim.lib.convolve_irf(pl, dtime)   

# Aligns sim either by max. (align_by_max) or steep (align_by_steep)
aligned_sims = kin_kit.align_by_max(sims, dtime, value = 0.5*ns)

fig = None
for i in range(len(aligned_sims)):
    color = ['#E24A33', '#348ABD'][i]
    fig = artists.plot3scales.plot(dtime, aligned_sims[i], system, t_dict=to, ivtype='none', 
                         annotate=True, fig = fig, color=color, 
                         linewidth=6, opaq=0.4,
                         mlabel='hetero sim')
