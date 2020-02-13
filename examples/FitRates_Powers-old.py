"""
Fit with Rates - Heterostructure

Example file for fitting multiple data traces of the same compound, but different
excitation powers. A set of rate constants is sought that will successfully
simulate the data at all excitation powers.

"""
import os
from TRPL import sim
from TRPL import fit
from TRPL import artists as art
from TRPL import kit as TRPL_kit
import matplotlib.pyplot as plt
from TRPL import data
from TRPL.units import *
import numpy as np
import scipy as sp
import time

doFit = False 
normalize = True #use false for power study, True for single heterostructure, irrelevant for single structure

#--- Creating Time Object
to = sim.time.linear(N=1000, subsample=1, pulse_window = 2*ps, numcycles=1000)
dtime = to['array'][::to['subsample']]

#--- Creating Data Object(s)
all_data = []; keys = []; powers = []
#dir_path = os.path.dirname(os.path.realpath(__file__))
dir_path = r'C:\Users\Natalia Spitha\wright_stuff_chem@wisc.edu\PL Microscope\Data\Natalia\n_series\data'
subfolder = r'Ea' # set subfolder='' if no subfolder
sourcefolder = os.path.join(dir_path, subfolder) # important: any saved data will go to this folder
destfolder = r'C:\Users\Natalia Spitha\wright_stuff_chem@wisc.edu\PL Microscope\Data\Natalia\n_series\Fits\Powers'


filenames = [
             'Ea_n=1_0.045_0_TRPL.asc',
#             'Ea_n=1_0.424_1_TRPL.asc',
#             'Ea_n=1_0.963_1_TRPL.asc',
             ]

powers = []
species_sets = []

#title = filenames[0].split('_')[1] + ' (' + 'sample ' + filenames[0].split('_')[0] + ')'
title = subfolder + '_3powers_log'

for f, filename in enumerate(filenames):
    key = filename.split('_')[2] # can be anything, really
    keys.append(key)
    p = os.path.join(sourcefolder, filename)
    all_data.append(data.lib.data_from_TCSPC(p, key = key, t_zero=0.582259, 
                                             skip_h=200+154, skip_f=884,
                                             metadata=True, weigh_by_coll=True))
    
    all_data[f].pulse_power = float(filename.split('_')[2])*1e-6
    powers.append(float(filename.split('_')[2])*1e-6)
    

#--- Interpolate data x and y on time array
all_y = []
for do in all_data: # do is a data object
    ydata = do.dark_subtraction(os.path.join(dir_path, 'dark.asc'))
    ydata = do.interp(dtime)
    ydata = do.remove_zeros(ydata)
    all_y.append(ydata)
if len(all_y)==1:
    all_y = all_y[0]
all_y = np.array(all_y)

#--- Parameters of simulation and initializing of system
initparams = {
        'k_ann': 9.77e8,
        'k_dis': 6.92e8,
        'k_trx': 0,
        'k_eth': 1e4,
        'k_rec': 2.14e5,
        'N_trp': 2e3,
        'cs': 0.19,
        } 

#--- Fitting boundaries
bounds = {
    'k_ann': (6e7, 6e9),
    'k_dis': (8e8, 2e10),
    'k_trx': (5e3, 1e7),
    'k_eth': (2e2, 3e4),
    'k_rec': (2e3, 3e5),
    'N_trp': (1e2, 5e4),
    'cs': (0.22, 0.48)
    }  

#--- Bounds as used in differential_evolution function
boundtuples = list(bounds.values())
#--- Bounds as used in least_squares function
minbounds, maxbounds = [i for i in zip(*boundtuples)]

#--- Create Excitation object and initial system
light = sim.lib.Excitation()
system = sim.systems.MonoHTrap(**initparams)

# Conditions for fitting algorithm
sim_comp_args = fit.lib.sac_args(
        varparamkeys=bounds.keys(),
        system = system, 
        data_arrays = all_y,
        to = to,
        light = light, 
        powers = powers,
        irf_fwhm = 40*ps, 
        roll_value=0, # in time units
        comparison='linear', # 'linear' or 'log'
        absolute=True,
        limits=None, # None or list of two time values
        norm = normalize,
        roll_criterion='max', #'max' or 'steep'
        maxavgnum=5,
        )


if doFit:
    time_start = time.clock()
    opt = sp.optimize.differential_evolution(fit.lib.simulate_and_compare,
                                             bounds = boundtuples, 
                                             args = (sim_comp_args))    
    fitparams = opt.x
    time_elapsed = time.clock() - time_start
    print('Fitting took %0.5f seconds'%(time_elapsed))
    fitparamdict = TRPL_kit.dict_from_list(fitparams, bounds.keys())
    TRPL_kit.printparamsexp(fitparamdict)

else: 
    fitparamdict = initparams

#--- Redefine system and simulation based on fitparams
system.update(**fitparamdict)
pulse = light.pulse_params()
for i, power in enumerate(powers):
            pulse['power'] = power
            light = sim.lib.Excitation(pulse=pulse)
            transient, converged = sim.lib.refined_simulation(system, to, light,
                                                      N_coarse=N_coarse)
            species_sets.append(transient)
            pl_at_this_power = system.PLsig(transient)
            
            if i == 0:
                pl = pl_at_this_power
            else:
                pl = np.vstack((pl, pl_at_this_power))
            
#pl = system.PLsig(transient)
sims = sim.lib.convolve_irf(pl, dtime)   

all_y = TRPL_kit.make_2d(all_y)
sims = TRPL_kit.make_2d(sims)
species_sets = np.array(species_sets)

#--- Re-align simulation and experimental traces for display and saving.
# can use align_by_steep or align_by_max
aligned_data = TRPL_kit.align_by_steep(all_y, dtime, value=0.5*ns)
aligned_sims = TRPL_kit.align_by_steep(sims, dtime, value = 0.5*ns)

colors=['#E24A33', '#348ABD', '#988ED5', '#8EBA42', '#FBC15E', '#777777', '#FFB5B8']

for i in range(len(aligned_data)):
    aligned_data[i] /= max(aligned_data[np.argmax(powers)])
    aligned_sims[i] /= max(aligned_sims[np.argmax(powers)])
    newFig=i==0
    art.plot3scales.plot(dtime, aligned_data[i], sys_obj=None, t_dict=None, 
                     annotate=False, newFig=newFig, linewidth=1, color= colors[i],
                     mlabel=keys[i]+' data', norm=normalize)
    
    art.plot3scales.plot(dtime, aligned_sims[i], system, t_dict=to, ivtype='none', 
                         annotate=True, newFig=False, linewidth=6, opaq=0.4, color=colors[i],
                         mlabel=keys[i] + ' sim', norm=normalize)

#--- Saving data into dictionary
trace_dict = {}
trace_dict['time (ns)'] = dtime/ns
for i, key in enumerate(keys):
    trace_dict[key+'_data'] = aligned_data[i]
    trace_dict[key+'_sims'] = aligned_sims[i]
    for i, species in enumerate(species_sets[i]):
        trace_dict[key + '_' + system.populations[i]] = species

saveparam_dict = system.params()

if doFit:
    saveparam_dict['opt.fun'] = opt.fun
else:
    saveparam_dict['opt.fun'] = 'fit not attempted'


#TRPL_kit.save_fit(title, trace_dict, saveparam_dict, destfolder) # comment out if saving not desired
#art.plot3scales.show_species(dtime, species_sets, system, False)
#plt.savefig(os.path.join(destfolder, title)+'_species.svg', format='svg')
#plt.savefig(os.path.join(destfolder, title)+'._species.png', format='png')
#art.plot3scales.show_species(dtime, species_sets, system, True)
#plt.savefig(os.path.join(destfolder, title)+'_speciesnorm.svg', format='svg')
#plt.savefig(os.path.join(destfolder, title)+'._speciesnorm.png', format='png') 
