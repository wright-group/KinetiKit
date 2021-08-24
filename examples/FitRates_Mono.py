"""
Fit with monolithic models (usually one power and PL output)

Example script for fitting TRPL data with a simulation of the Mono type, 
plotting, and saving it.

See QUICK NOTES in README file () (github.com/wright-group/KinetiKit) for help
"""

import os
import time

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from KinetiKit import sim, fit, data, artists
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
key = 'perovskite' # identifier for your data (used in plot legend)
dir_path = os.path.dirname(os.path.realpath(__file__)) # directory path
subfolder = 'ex_data' # enter '' if no subfolder
filename = r'sample_TRPL.asc'

#--- file details:
skip_header = 154; skip_footer = 884 # see data file for how many lines to skip
collection_time = 900 # data collection time in seconds

#--- subtract dark counts? (ignore all but next line if False)
sub_dark_counts = True
dark_path = os.path.join(dir_path, subfolder, 'dark.asc')
dark_skip_header = 354; dark_skip_footer = 884
dark_collection_time = 300 

#--- output 
# CALL save_all() ON THE CONSOLE AFTER FITTING IS COMPLETE,
# TO SAVE PARAMETERS AND PLOTS
destpath = os.path.join(dir_path, 'ex_output', 'FitRates_Mono')
output_name = key # can be anything, DO NOT include extension


"""
Define Time and Excitation Parameters
"""
#--- arguments of sim.lib.Excitation()
pulse_power = 1000*nW
reprate = 80 * MHz
pulse_fwhm = 100*fs
pulse_wavelength = 400*nm
#--- arguments of sim.time.linear()
N= 1000
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
irf_args = {'irf_type': 'GaussDiff',
            'weighted' : True,
            'fwhm': 55 * ps,
            'tau': 800 *ps,
            'b': 0.07,
            'tau_wt': 10 *ps
            }

"""
Choose System type, Initial parameters, and search boundaries. 
See QUICK NOTES > Note on Setting and Fitting parameters in README file
"""
system = sim.systems.Mono()

initparams = {
    'k_ann': 1.25e8,
    'k_dis': 2.01e9,
    'k_rec': 1.3e3,
    'cs': 0.5,
    }

bounds = {
        'k_ann': (1e6, 1e10),
        'k_dis': (5e6, 3.2e10), 
        'k_rec': (1e3, 1e8),
        } 

"""
Fitting preferences
"""
doFit = True # if False, system will be modeled with initparams
doLS = True # whether to refine the optimization via a local least-squares 
            # fitting (and obtain error estimates). Ignore if doFit = False
settings['display_counter'] = True # display counter showing search iteration

#--- arguments of sim.fit.simulate_and_compare() -- see docstring
comparison_type = 'linear' # "linear" of "log" comparison betw. data and sim.
absolute = True # return avg. of absolute vs. relative differences
norm = True # False recommended for simultaneous multi-power fitting;
            # see docstring
roll_criterion = 'max' # 'max' or 'steep': : whether to align data and
                       # simulation based on their maximum or steepest point
avgnum = 5 # how many points to consider to determine the max/steepest point

###############################################################################
# End of editable section
###############################################################################

#--- Creating Time Object
to = sim.time.linear(N=N, period=period, subsample = subsample)
dtime = to['array'][::to['subsample']]

#--- Creating Data Object(s)
file_path = os.path.join(dir_path, subfolder,  filename)
do = data.lib.data_from_SPCM(file_path, key = key,
                                    skip_h=skip_header, skip_f=skip_footer,
                                    weigh_by_coll=True, coll=collection_time)

if sub_dark_counts:
    dark_counts = data.lib.data_from_SPCM(dark_path,
                                    skip_h=dark_skip_header, 
                                    skip_f=dark_skip_footer,
                                    weigh_by_coll=True,
                                    coll = dark_collection_time)
    do.dark_subtract(dark_counts, method='elementwise')

# Data is interpolated to match same time axis as simulation
do.interp(dtime)

all_y = do.y

#--- Parameters of simulation and initializing of system

system.update(**initparams)
boundtuples = list(bounds.values()) # Bounds as used in differential_evolution function

#--- Create Excitation object
light = sim.lib.Excitation(pulse={'power':pulse_power,
                                  'reprate': reprate,
                                  'fwhm' : pulse_fwhm,
                                  'wavelength': pulse_wavelength})

# Conditions for fitting algorithm
conditions = fit.lib.sac_args(
        varparamkeys=bounds.keys(),
        system = system, 
        data_arrays = all_y,
        to = to,
        light = light,
        irf_args = irf_args,
        N_coarse = N_coarse,  
        comparison=comparison_type, # 'linear' or 'log'
        absolute=absolute,
        limits=limits, # None or list of two time values
        norm=norm,
        roll_criterion=roll_criterion, 
        maxavgnum=avgnum
        )

if doFit:
    time_start = time.time()
    counter = 0
    # First perform a global search using Differential Evolution
    opt_DE = sp.optimize.differential_evolution(fit.lib.simulate_and_compare,
                                              bounds= boundtuples, 
                                              args= conditions,
                                              )
    if doLS:
        # Fine-tune with a least-squares fit to determine curvature 
        # of parameter space
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
    
    time_elapsed = time.time() - time_start
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
transient, converged = sim.lib.refined_simulation(system, to, light, N_coarse = N_coarse) # system re-simulated with fit parameters
species_set = transient # set of population's time evolution
pl = system.PLsig(transient) # PLsig function applied to population arrays
sims = sim.lib.convolve_irf(pl, dtime, irf_args)    # PL signal convolved with IRF

# alignment
if roll_criterion == 'max':
    aligned_data = kin_kit.align_by_max(all_y, dtime, value = align_to, avgnum=avgnum)
    aligned_sims = kin_kit.align_by_max(sims, dtime, value = align_to, avgnum=avgnum)
elif roll_criterion == 'steep':
    aligned_data = kin_kit.align_by_steep(all_y, dtime, value = align_to, avgnum=avgnum)
    aligned_sims = kin_kit.align_by_steep(sims, dtime, value = align_to, avgnum=avgnum)


#--- Plot
fig = artists.plot3scales.plot(dtime, aligned_data, sys_obj=None, t_dict=None,
                     unit = time_unit, annotate=False, fig=None, linewidth=1, 
                     mlabel=key+' data')
    
artists.plot3scales.plot(dtime, aligned_sims, system, t_dict=to, ivtype='none', 
                     unit = time_unit, annotate=True, fig=fig, ResetColorCyc=True, 
                     linewidth=6, opaq=0.4,
                     mlabel=key+' sim')

#--- Saving data and parameters into dictionary
trace_dict = {}
trace_dict['time (%s)'%time_unit] = dtime/units[time_unit]

trace_dict[key+'_data'] = aligned_data
trace_dict[key+'_sims'] = aligned_sims

for i, species in enumerate(species_set):
    trace_dict[key + '_' + system.populations[i]] = species
        
saveparam_dict = kin_kit.saveparam_dict(system, opt_DE, opt_LS, bounds, conditions, doFit)

def save_all(): 
    # save all images, parameters, and traces, in image and CSV files
    kin_kit.save_fit(output_name, trace_dict, saveparam_dict, destpath)
    plot_species(save=True)
    
def plot_species(save=False):
    artists.plot3scales.show_species(dtime, [species_set], system, False, 
                                     save=save, filename=output_name, destfolder=destpath)
