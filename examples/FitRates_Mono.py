"""
Fit with monolithic models (usually one power and PL output)

Example file for fitting a data object with a simulation of the Mono type. 
Requires numpy, scipy, matplotlib, and time packages along with KinetiKit package
"""

import os
import time

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from KinetiKit import sim, fit, data, artists
from KinetiKit import kit as kin_kit
from KinetiKit.units import nW, uW, MHz, nm, ps, ns
from KinetiKit.settings import settings


doFit = True
doLS = True # whether to refine the optimization via a local least-squares fitting (and obtain error estimates)

settings['display_counter'] = True # display counter showing search iteration

#--- Creating Time Object
to = sim.time.linear(N=1000, period=12.5*ns)
dtime = to['array'][::to['subsample']]

#--- Create System Instance (select model from sim.systems)
system = sim.systems.Mono()

#--- Creating Data Object(s)
dir_path = os.path.dirname(os.path.realpath(__file__)) # this will read files from the path this file is in
#dir_path = r'C:\...' # alternatively you can specify the path from which to read files
subfolder = r'ex_data' # set subfolder='' if no subfolder
destfolder = os.path.join(dir_path, r'ex_output', r'FitRates_Mono') # folder where any output files will be saved

filename = r'pure_metadata_TRPL.asc'

key = 'perovskite' # this will be the legend label, can be anything
power = 1000*nW #Excitation power

dark_path = os.path.join(dir_path, subfolder, 'dark.asc')
dark_counts = data.lib.data_from_SPCM(dark_path,
                                    skip_h=353, skip_f=884,
                                    metadata=True, weigh_by_coll=True)

file_path = os.path.join(dir_path, subfolder,  filename)
do = data.lib.data_from_SPCM(file_path, key = key,
                                    skip_h=353, skip_f=884,
                                    weigh_by_coll=True, metadata=True, coll=60)


# Various transformations performed on data
do.dark_subtract(dark_counts, method='elementwise')
do.interp(dtime)
do.remove_zeros()

all_y = do.y

#--- Parameters of simulation and initializing of system
# Note that, as long as bounds have been set for a parameter, its value given
# in initparams does NOT affect the Differential Evolution search
initparams = {
    'k_ann': 1.25e8,
    'k_dis': 2.01e9,
    'k_rec': 1.3e3,
    'cs': 0.5,
    }

system.update(**initparams)

#--- Fitting boundaries
bounds = {
        'k_ann': (1e6, 1e10),
        'k_dis': (5e6, 3.2e10), 
        'k_rec': (1e3, 1e8),
        }  


boundtuples = list(bounds.values()) # Bounds as used in differential_evolution function

#--- Create Excitation object
light = sim.lib.Excitation(pulse={'power':power,
                                  'reprate': 80*MHz,
                                  'wavelength': 400*nm})

# Conditions for fitting algorithm
conditions = fit.lib.sac_args(
        varparamkeys=bounds.keys(),
        system = system, 
        data_arrays = all_y,
        to = to,
        light = light, # this and arguments above it do not need to be changed
        irf_fwhm = 40*ps, # Instrument Response Function
        N_coarse = 500,  
        comparison='linear', # 'linear' or 'log'
        absolute=True,
        limits=None, # None or list of two time values
        norm=True,
        roll_criterion='max', #'max' or 'steep': whether to align data and 
                              # simulation based on their maximum or steepest point
        maxavgnum=5
        )

if doFit:
    time_start = time.clock()
    counter = 0
    # First perform a global search using Differential Evolution
    opt_DE = sp.optimize.differential_evolution(fit.lib.simulate_and_compare,
                                              bounds= boundtuples, 
                                              args= conditions,
                                              )
    if doLS:
        # Fine-tune with a least-squares fit to determine curvature of parameter space
        counter = 0
        opt_LS = fit.lib.fit_leastsq(fit.lib.simulate_and_compare, 
                                        p0 = opt_DE.x, 
                                        args=tuple(list(conditions[:-1]) + [False]))
        fitparams = opt_LS[0]
        errordict = kin_kit.dict_from_list(opt_LS[1], bounds.keys())
    else:
        opt_LS = None
        errordict = None
        fitparams = opt_DE.x
    
    time_elapsed = time.clock() - time_start
    print('Fitting took %0.5f seconds'%(time_elapsed))
    fitparamdict = kin_kit.dict_from_list(fitparams, bounds.keys())
    
else: 
    fitparamdict = initparams
    errordict=None
    opt_DE=None
    opt_LS=None
    
#--- Redefine system and simulation based on fitparams
kin_kit.printparamsexp(fitparamdict, errordict) # display fit parameters

system.update(**fitparamdict) # system updated with fit parameters
transient, converged = sim.lib.refined_simulation(system, to, light) # system re-simulated with fit parameters
species_set = transient # set of population's time evolution
pl = system.PLsig(transient) # PLsig function applied to population arrays
sims = sim.lib.convolve_irf(pl, dtime)   # PL signal convolved with IRF

# Aligns data with sim either by max. (align_by_max) or steep (align_by_steep)
align_to = 0.5*ns #value at which data and simulation are aligned for saving/diplay
aligned_data = kin_kit.align_by_max(all_y, dtime, value = align_to, avgnum=5)
aligned_sims = kin_kit.align_by_max(sims, dtime, value = align_to, avgnum=5)

#--- Plot
fig = artists.plot3scales.plot(dtime, aligned_data, sys_obj=None, t_dict=None,
                     annotate=False, fig=None, linewidth=1, 
                     mlabel=key+' data')
    
artists.plot3scales.plot(dtime, aligned_sims, system, t_dict=to, ivtype='none', 
                     annotate=True, fig=fig, ResetColorCyc=True, 
                     linewidth=6, opaq=0.4,
                     mlabel=key+' sim')

#--- Saving data and parameters into dictionary
trace_dict = {}
trace_dict['time (ns)'] = dtime/ns

trace_dict[key+'_data'] = aligned_data
trace_dict[key+'_sims'] = aligned_sims

for i, species in enumerate(species_set):
    trace_dict[key + '_' + system.populations[i]] = species
        
saveparam_dict = kin_kit.saveparam_dict(system, opt_DE, opt_LS, bounds, conditions, doFit)

def save_all(): 
    # save all images, parameters, and traces, in image and CSV files
    kin_kit.save_fit(key, trace_dict, saveparam_dict, destfolder)
    plot_species(save=True)
    
def plot_species(save=False):
    artists.plot3scales.show_species(dtime, [species_set], system, False, 
                                     save=save, filename=key, destfolder=destfolder)
