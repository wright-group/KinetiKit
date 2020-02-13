"""
Simultaneous fit example with same model at two different powers

See QUICK NOTES in README file () (github.com/wright-group/KinetiKit) for help
note: if you run doFit on this program as originally in the repository, 
fitting may take up to 10 minutes

"""

import os
import time

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from KinetiKit import sim, fit, data
from KinetiKit import artists as art
from KinetiKit import kit as kin_kit
from KinetiKit.units import units, nW, uW, MHz, nm, ps, ns, fs
from KinetiKit.settings import settings


###############################################################################
# Start of editable section
###############################################################################

"""
Define relevant filenames. See QUICK NOTES > Note on File paths 
"""
#--- data to be imported 
keys = ['45 nW', '400 nW'] # identifiers for your data (used in plot legend)
dir_path = os.path.dirname(os.path.realpath(__file__)) # directory path
subfolders = [r'ex_data/low', r'ex_data/mid'] # enter '' if no subfolder
filenames = [r'low_power.asc', 
             r'mid_power.asc']

#--- file details:
skip_header = 354; skip_footer = 884 # see data file for how many lines to skip
collection_times = [600, 180] # data collection time in seconds for each file

#--- subtract dark counts? (ignore all but next line if False)
sub_dark_counts = True
dark_path = os.path.join(dir_path, 'ex_data', 'dark.asc')
dark_skip_header = 354; dark_skip_footer = 884
dark_collection_time = 300 

#--- output 
# CALL save_all() ON THE CONSOLE AFTER FITTING IS COMPLETE,
# TO SAVE PARAMETERS AND PLOTS
destpath = os.path.join(dir_path, 'ex_output', 'FitRates_Multipower')
output_name = 'multipower' # can be anything, DO NOT include extension


"""
Define Time and Excitation Parameters
"""
#--- arguments of sim.lib.Excitation()
pulse_powers = [45*nW, 400*nW]
reprate = 80 * MHz
pulse_fwhm = 100*fs
pulse_wavelength = 400*nm
#--- arguments of sim.time.linear()
N = 1000
time_unit = 'ns'
period = 1/reprate # simulation time axis should span (1/reprate) of laser
subsample = 1 
#--- other
limits = None # set to [start*ns, end*ns] when data time < 1/reprate
N_coarse = 500 # number of time points for coarse simulation
align_to = 0.5*ns # value to which data and simulation are aligned for saving

"""
Define Instrument Response Function FWHM (Note: currently, only simulated
FWHM are supported)
"""
irf_fwhm = 40*ps

"""
Choose System type, Initial parameters, and search boundaries. 
See QUICK NOTES > Note on Setting and Fitting parameters in README file
"""
system = sim.systems.MonoRecX()

initparams = {
        'k_ann': 1.0e9,#1001463356.9975033,
        'k_dis': 9e9,#9990872345.442875,
        'k_rec': 5.299e7,#52987975.39923163,
        'cs': 0.5,
    }  

bounds = {
        'k_ann': (1e7, 1e10),
        'k_dis': (1e9, 9.9e9), 
        'k_rec': (1e5, 1e8),
        } 

"""
Fitting preferences
"""
doFit = True # if False, system will be modeled with initparams
doLS = True # whether to refine the optimization via a local least-squares 
            # fitting (and obtain error estimates). Ignore if doFit = False
settings['display_counter'] = True # display counter showing search iteration

#--- arguments of sim.fit.simulate_and_compare() -- see docstring
comparison_type = 'log' # "linear" of "log" comparison betw. data and sim.
absolute = True # return avg. of absolute vs. relative differences
norm = False # False recommended for simultaneous multi-power fitting;
            # see docstring
roll_criterion = 'steep' # 'max' or 'steep': : whether to align data and
                       # simulation based on their maximum or steepest point
avgnum = 5 # how many points to consider to determine the max/steepest point

###############################################################################
# End of editable section
###############################################################################

#--- Creating Time Object
to = sim.time.linear(N=N, period=period, subsample = subsample)
dtime = to['array'][::to['subsample']]


#--- Creating Data Object(s)
all_y = []
for i, filename in enumerate(filenames):
    file_path = os.path.join(dir_path, subfolders[i],  filename)
    do = data.lib.data_from_SPCM(file_path, key = keys[i],
                                        skip_h=skip_header, skip_f=skip_footer,
                                        weigh_by_coll=True, coll=collection_times[i])
    if sub_dark_counts:
        dark_counts = data.lib.data_from_SPCM(dark_path,
                                        skip_h=dark_skip_header, 
                                        skip_f=dark_skip_footer,
                                        weigh_by_coll=True,
                                        coll = dark_collection_time)
        do.dark_subtract(dark_counts, method='average')
    
    # Data is interpolated to match same time axis as simulation
    do.interp(dtime)
    
    all_y.append(do.y)
all_y = np.array(all_y)

#--- Parameters of simulation and initializing of system

system.update(**initparams)
boundtuples = list(bounds.values()) # Bounds as used in differential_evolution function

#--- Create Excitation object
light = sim.lib.Excitation(pulse={'reprate': reprate,
                                  'fwhm' : pulse_fwhm,
                                  'wavelength': pulse_wavelength})

# Conditions for fitting algorithm
conditions = fit.lib.sac_args(
        varparamkeys=bounds.keys(),
        system = system, 
        data_arrays = all_y,
        to = to,
        light = light,
        powers = pulse_powers, # each power in this list will update the light
                             # object appropriately
        irf_fwhm = irf_fwhm,
        N_coarse = N_coarse,  
        comparison=comparison_type, # 'linear' or 'log'
        absolute=absolute,
        limits=limits, # None or list of two time values
        norm=norm,
        roll_criterion=roll_criterion, 
        maxavgnum=avgnum
        )

if doFit:
    time_start = time.clock()
    counter = 0
    # First perform a global search using Differential Evolution
    opt_DE = sp.optimize.differential_evolution(fit.lib.simulate_and_compare,
                                              bounds= boundtuples, 
                                              args= conditions,
                                              )
    print("differential evolution search complete")
    if doLS:
        # Fine-tune with a least-squares fit to determine curvature 
        # of parameter space
        counter = 0
        opt_LS = fit.lib.fit_leastsq(fit.lib.simulate_and_compare, 
                                        p0 = opt_DE.x, 
                                        args=tuple(list(conditions[:-1]) + [False]))
        fitparams = opt_LS[0]
        errordict = kin_kit.dict_from_list(opt_LS[1], bounds.keys())
        print("least squares search complete")
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
species_sets = []
for i, power in enumerate(pulse_powers):
    light = light.updated_with(pulse={'power':power})
    # print(i, 'pulse', light.pulse, 'cw', light.cw, 'numc', light.numcycles)
    # print('system params', system.params(), 'time', to.items())
    transient, converged = sim.lib.refined_simulation(system, to, light,
                                              N_coarse=N_coarse)
    # print(converged)
    species_sets.append(transient)
    pl_at_this_power = system.PLsig(transient)
    
    if i == 0:
        pl = pl_at_this_power
    else:
        pl = np.vstack((pl, pl_at_this_power))
                
sims = sim.lib.convolve_irf(pl, dtime)   # PL signal convolved with IRF
# all_y = kin_kit.make_2d(all_y)
# sims = kin_kit.make_2d(sims)
species_sets = np.array(species_sets)

# alignment
if roll_criterion == 'max':
    aligned_data = kin_kit.align_by_max(all_y, dtime, value = align_to, avgnum=avgnum)
    aligned_sims = kin_kit.align_by_max(sims, dtime, value = align_to, avgnum=avgnum)
elif roll_criterion == 'steep':
    aligned_data = kin_kit.align_by_steep(all_y, dtime, value = align_to, avgnum=avgnum)
    aligned_sims = kin_kit.align_by_steep(sims, dtime, value = align_to, avgnum=avgnum)


#--- Plot
colors = art.plotparams.colors
fig = None
for i in range(len(aligned_data)):
    aligned_data[i] /= max(aligned_data[np.argmax(pulse_powers)])
    aligned_sims[i] /= max(aligned_sims[np.argmax(pulse_powers)])
    fig = art.plot3scales.plot(dtime, aligned_data[i], sys_obj=None, t_dict=None, 
                     annotate=False, fig = fig, linewidth=1, color= colors[i],
                     mlabel=keys[i]+' data', norm=norm)
    
    art.plot3scales.plot(dtime, aligned_sims[i], system, t_dict=to, ivtype='none', 
                         annotate=True, fig=fig, linewidth=6, opaq=0.4, color=colors[i],
                         mlabel=keys[i] + ' sim', norm=norm)


#--- Saving data and parameters into dictionary
trace_dict = {}
trace_dict['time (%s)'%time_unit] = dtime/units[time_unit]

for i, key in enumerate(keys):
    trace_dict[key+'_data'] = aligned_data[i]
    trace_dict[key+'_sims'] = aligned_sims[i]
    for i, species in enumerate(species_sets[i]):
        trace_dict[key + '_' + system.populations[i]] = species
        
saveparam_dict = kin_kit.saveparam_dict(system, opt_DE, opt_LS, bounds, conditions, doFit)

def save_all(): 
    # save all images, parameters, and traces, in image and CSV files
    kin_kit.save_fit(output_name, trace_dict, saveparam_dict, destpath)
    plot_species(save=True)
    
def plot_species(save=False):
    artists.plot3scales.show_species(dtime, species_sets, system, False, 
                                     save=save, filename=key, destfolder=destpath)
