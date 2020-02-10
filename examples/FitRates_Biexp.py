"""
Fit with Rates_Mono model

Example file for fitting a data object with a biexponential function. 
Requires numpy, scipy, matplotlib, and time packages along with TRPL package
"""


import os
import time

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from KinetiKit import sim, fit, data
from KinetiKit import artists as art
from KinetiKit import kit as kin_kit
from KinetiKit.units import ns, ps


doFit = True

#--- Creating Time Object
to = sim.time.linear(N=1000); dtime = to['array'][::to['subsample']]

#--- Create system instance
system = sim.systems.Biexp()

#--- Creating Data Object(s)
dir_path = os.path.dirname(os.path.realpath(__file__))

filename = 'pure_metadata_TRPL.asc'
    
subfolder = 'ex_data'
destfolder = os.path.join(dir_path, r'ex_output')

all_data = []; keys = []; powers = []
#filenames = [r'Xb_n=0_0.045_0_TRPL.asc']

key = 'FitRates_Biexp'
keys.append(key)
p = os.path.join(dir_path, subfolder,  filename)
do = data.lib.data_from_SPCM(p, key = key, skip_h=200+154, skip_f=884,
                                         weigh_by_coll=True, metadata=True)
    
#--- Transformations on data
do.dark_subtract(os.path.join(dir_path, subfolder, 'dark.asc'), 300)
do.interp(dtime)
do.remove_zeros()
all_y = do.y

#--- Parameters of simulation and initializing of system

initparams = {
           'A1': 1.9974e-1,
           'tau1': 1.774e-9,
           'tau2': 1e-10,
                }

system.update(**initparams)

#--- Fitting boundaries
bounds = {
        'A1': (0, 1),
        'tau1': (1e-10, 1e-6),
        'tau2': (1e-10, 1e-6),
        }    

#--- Bounds as used in differential_evolution function
boundtuples = list(bounds.values())
#--- Bounds as used in least_squares function
minbounds, maxbounds = [i for i in zip(*boundtuples)]

# Conditions for fitting algorithm
sim_comp_args = fit.lib.sac_args(
        varparamkeys=bounds.keys(),
        system = system, 
        data_arrays = all_y,
        to = to,
        light = None,
        powers=None,
        irf_fwhm = 40*ps, 
        roll_value=0, # in time units
        comparison='linear', # 'linear' or 'log'
        absolute=True,
        limits=None, # None or list of two time values
        norm=True,
        roll_criterion='steep', #'max' or 'steep'
        maxavgnum=5
        )

if doFit:
    time_start = time.clock()
    opt = sp.optimize.differential_evolution(fit.lib.simulate_and_compare,
                                             bounds= boundtuples, 
                                             args= sim_comp_args)
    
    time_elapsed = time.clock() - time_start
    fitparams = opt.x
    print('Fitting took %0.5f seconds'%(time_elapsed))
    fitparamdict = kin_kit.dict_from_list(fitparams, bounds.keys())
    
else: 
    fitparamdict = initparams
    
#--- Redefine system and simulation based on fitparams
kin_kit.printparamsexp(fitparamdict) # display fit parameters

system.update(**fitparamdict)
pl, converged = sim.lib.simulate_func(system, dtime)
sims = sim.lib.convolve_irf(pl, dtime)   

# Aligns data with sim either by max. or steep
aligned_data = kin_kit.align_by_steep(all_y, dtime, value = 0.5*ns)
aligned_sims = kin_kit.align_by_steep(sims, dtime, value = 0.5*ns)

#--- Plot
fig=None
fig = art.plot3scales.plot(dtime, aligned_data, sys_obj=None, t_dict=None,
                     annotate=False, fig=fig, linewidth=1, 
                     #color='#348ABD',
                     mlabel=key+' data')
    
art.plot3scales.plot(dtime, aligned_sims, system, t_dict=to, ivtype='none', 
                     annotate=True, fig=fig, ResetColorCyc=True, 
                     #color = '#348ABD',
                     linewidth=6, opaq=0.4,
                     mlabel=key+' sim')

#--- Saving data and parameters into dictionary
trace_dict = {}
trace_dict['time (ns)'] = dtime/ns
for i, key in enumerate(keys):
    trace_dict[key+'_data'] = aligned_data
    trace_dict[key+'_sims'] = aligned_sims

if doFit:            
    saveparam_dict = kin_kit.saveparam_dict(system, opt, bounds, sim_comp_args, doFit)
    
def save():
    if not os.path.exists(destfolder):
        os.makedirs(destfolder)
        print("Created %s directory"%destfolder)
        
    kin_kit.save_fit(key, trace_dict, saveparam_dict, destfolder)
