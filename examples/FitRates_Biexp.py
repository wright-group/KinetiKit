"""
Fit with Rates_Mono model

Example file for fitting a data object with a biexponential function. 
Requires numpy, scipy, matplotlib, and time packages along with TRPL package
"""


import os
import time

import scipy as sp

from KinetiKit import sim, fit, data
from KinetiKit import artists as art
from KinetiKit import kit as kin_kit
from KinetiKit.units import ns, ps
from KinetiKit.settings import settings

#--- data to be imported 
key = r'FitRates_Biexp' # identifiers for your data (used in plot legend)
dir_path = os.path.dirname(os.path.realpath(__file__)) # directory path
subfolder = r'ex_data' # enter '' if no subfolder
filename = r'sample_TRPL.asc'

#--- file details:
skip_header = 154; skip_footer = 884 # see data file for how many lines to skip

#--- subtract dark counts? (ignore all but next line if False)
sub_dark_counts = True
dark_path = os.path.join(dir_path, 'ex_data', 'dark.asc')
dark_skip_header = 354; dark_skip_footer = 884
dark_collection_time = 300 

#--- output 
# CALL save_all() ON THE CONSOLE AFTER FITTING IS COMPLETE,
# TO SAVE PARAMETERS AND PLOTS
destpath = os.path.join(dir_path, 'ex_output/FitRates_Biexp')
output_name = 'biexponential' # can be anything, DO NOT include extension

"""
Define Time and Excitation Parameters
"""
#--- arguments of sim.time.linear()
N = 1000
time_unit = 'ns'
period = 12.5*ns # simulation time axis should span (1/reprate) of laser
subsample = 1 
N_coarse = 500

#--- other
limits = None # set to [start*ns, end*ns] when data time < 1/reprate
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
system = sim.systems.Biexp()
initparams = {
           'A1': 1.9974e-1,
           'tau1': 1.774e-9,
           'tau2': 1e-10,
           'offset': 0
                } 
bounds = {
        'A1': (0, 1),
        'tau1': (1e-10, 1e-6),
        'tau2': (1e-10, 1e-6),
        'offset': (0,1)
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
norm = True # False recommended for simultaneous multi-power fitting;
            # see docstring
roll_criterion = 'steep' # 'max' or 'steep': : whether to align data and
                       # simulation based on their maximum or steepest point
avgnum = 5 # how many points to consider to determine the max/steepest point

###############################################################################
# End of editable section
###############################################################################

#--- Creating Time Object
to = sim.time.linear(N=N, period=period, subsample=subsample)
dtime = to['array'][::to['subsample']]

#--- Creating Data Object(s)
p = os.path.join(dir_path, subfolder,  filename)
do = data.lib.data_from_SPCM(p, key = key, skip_h=skip_header, skip_f=skip_footer,
                                         weigh_by_coll=True)
    
#--- Transformations on data
if sub_dark_counts:
        dark_counts = data.lib.data_from_SPCM(dark_path,
                                        skip_h=dark_skip_header, 
                                        skip_f=dark_skip_footer,
                                        weigh_by_coll=True,
                                        coll = dark_collection_time)
        do.dark_subtract(dark_counts, method='average')
do.interp(dtime)
all_y = do.y

#--- Updating system with initial parameters
system.update(**initparams)

#--- Bounds as used in differential_evolution function
boundtuples = list(bounds.values())

# Conditions for fitting algorithm
conditions = fit.lib.sac_args(
        varparamkeys=bounds.keys(),
        system = system, 
        data_arrays = all_y,
        light=None,
        to = to,
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
kin_kit.printparamsexp(fitparamdict) # display fit parameters

system.update(**fitparamdict)
pl, converged = sim.lib.simulate_func(system, dtime)
sims = sim.lib.convolve_irf(pl, dtime)   

# Aligns data with sim either by max. or steep
if roll_criterion == 'max':
    aligned_data = kin_kit.align_by_max(all_y, dtime, value = align_to, avgnum=avgnum)
    aligned_sims = kin_kit.align_by_max(sims, dtime, value = align_to, avgnum=avgnum)
elif roll_criterion == 'steep':
    aligned_data = kin_kit.align_by_steep(all_y, dtime, value = align_to, avgnum=avgnum)
    aligned_sims = kin_kit.align_by_steep(sims, dtime, value = align_to, avgnum=avgnum)


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
trace_dict[key+'_data'] = aligned_data
trace_dict[key+'_sims'] = aligned_sims

if doFit:            
    saveparam_dict = kin_kit.saveparam_dict(system, opt_DE, opt_LS, bounds, conditions, doFit)
    
def save_all(): 
    # save all images, parameters, and traces, in image and CSV files
    kin_kit.save_fit(output_name, trace_dict, saveparam_dict, destpath)

