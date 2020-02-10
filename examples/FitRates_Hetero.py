"""
Fit with Rates - Heterostructure

Example file for fitting multiple data traces with a simulation of Hetero type, 
e.g. in a heterostructure with different emitting colors. 
Requires numpy, scipy, matplotlib, and time packages along with TRPL package

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

#note: if you run doFit on this program as originally in the repository, fitting may take up to 10 minutes
doFit = True # if false, a plot will be made using the defined initparams
doLS = False # whether to refine the optimization via a local least-squares fitting (and obtain error estimates)

settings['display_counter'] = True # display counter showing search iteration

#--- Creating Time Object
to = sim.time.linear(N=1000, period = 12.5*ns)
dtime = to['array'][::to['subsample']]

#--- Create System Instance (select model from sim.systems)
system = sim.systems.Hetero()

#--- Creating Data Object(s)

dir_path = os.path.dirname(os.path.realpath(__file__))
#dir_path = r'C:\Users\...'
subfolder = r'ex_data' # set subfolder='' if no subfolder
destfolder = os.path.join(dir_path, r'ex_output') # important: any saved data will go to this folder

all_y = []
commonkey = 'heterostructure'
keys = ['n=1', 'n=2']
power = 1000*nW
filenames = ['het_n=1_TRPL.asc', 
             'het_n=2_TRPL.asc']


dark_path = os.path.join(dir_path, subfolder, 'dark.asc')
dark_counts = data.lib.data_from_SPCM(dark_path,
                                    skip_h=353, skip_f=884,
                                    metadata=True, weigh_by_coll=True)

for f, filename in enumerate(filenames):
    key = keys[f]
    file_path = os.path.join(dir_path, subfolder, 'hetero', filename)
    do = data.lib.data_from_SPCM(file_path, key = key,
                                             skip_h=153, skip_f=884,
                                             coll=60)
    do.dark_subtract(dark_counts, method='elementwise')
    do.interp(dtime)
    do.remove_zeros()
    all_y.append(do.y)
    
all_y = np.array(all_y)

#--- Parameters of simulation and initializing of system
# Note that, as long as bounds have been set for a parameter, its value given
# in initparams does NOT affect the Differential Evolution search

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

system.update(**initparams)
#--- Fitting boundaries
bounds = {
        'k1_ann': (1e8, 2e10),
        'k1_dis': (1e8, 2e10),
        'k1_rec': (1e3, 1e7),
        'k2_ann': (1e8, 9e9),
        'k2_dis': (1e8, 2e10),
        'k2_rec': (1e3, 1e7),
        'k_xtr' : (1e7, 1e10),
}  

#--- Bounds as used in differential_evolution function
boundtuples = list(bounds.values())

#--- Create Excitation object and initial system
light = sim.lib.Excitation(pulse={'power':power,
                                  'reprate': 80*MHz,
                                  'wavelength': 400*nm})

# Conditions for fitting algorithm
conditions = fit.lib.sac_args(
        varparamkeys=bounds.keys(),
        system = system, 
        data_arrays = all_y,
        to = to,
        light = light, 
        irf_fwhm = 40*ps, 
        N_coarse = 500,
        comparison='log', # 'linear' or 'log'
        absolute=True,
        limits=None, # None or list of two time values
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
    
        fit.lib.fit_leastsq(fit.lib.simulate_and_compare, 
                                    p0 = opt_DE.x, 
                                    args=tuple(list(conditions[:-1]) + [False]))
        errordict = kin_kit.dict_from_list(opt_LS[1], bounds.keys())
        fitparams = opt_LS[0]
    
    else:
        opt_LS=None
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

system.update(**fitparamdict)
transient, converged = sim.lib.refined_simulation(system, to, light, N_coarse=conditions[7])
species_set = transient
pl = system.PLsig(transient)
sims = sim.lib.convolve_irf(pl, dtime)   

# Aligns data with sim either by max. (align_by_max) or steep (align_by_steep)
aligned_data = kin_kit.align_by_max(all_y, dtime, value=0.5*ns)
aligned_sims = kin_kit.align_by_max(sims, dtime, value = 0.5*ns)

fig = None
for i in range(len(aligned_data)):
    color = ['#E24A33', '#348ABD'][i]
    fig = artists.plot3scales.plot(dtime, aligned_data[i], sys_obj=None, t_dict=None, 
                     annotate=False, fig = fig, linewidth=1, color=color,
                     mlabel=keys[i]+' data')
    artists.plot3scales.plot(dtime, aligned_sims[i], system, t_dict=to, ivtype='none', 
                         annotate=True, fig = fig, color=color, 
                         linewidth=6, opaq=0.4,
                         mlabel=keys[i] + ' sim')

#--- Saving data into dictionary
trace_dict = {}
trace_dict['time (ns)'] = dtime/ns
for i, key in enumerate(keys):
    trace_dict[key+'_data'] = aligned_data[i]
    trace_dict[key+'_sims'] = aligned_sims[i]

saveparam_dict = kin_kit.saveparam_dict(system, opt_DE, opt_LS, bounds, conditions, doFit)

def save_all(): 
    # save all images, parameters, and traces, in image and CSV files
    kin_kit.save_fit(commonkey, trace_dict, saveparam_dict, destfolder)
    plot_species(save=True)
    
def plot_species(save=False):
    artists.plot3scales.show_species(dtime, [species_set], system, False, 
                                     save=save, filename=commonkey, destfolder=destfolder)