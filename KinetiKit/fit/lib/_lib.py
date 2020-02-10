"""
Fitting library 

Library of functions and classes that relate and compare data with simulations
for the purposes of fitting.

"""
import numpy as np
import scipy as sp

from KinetiKit import sim
from KinetiKit.units import units, ps
from KinetiKit.settings import settings, counter
import KinetiKit.kit as kin_kit


def elementwise_diff(data_arrays, sim_arrays, norm=True, comparison='linear',
                     absolute=False,):
    """
    Returns a list of arrays of the relative differences between pairs of 
    experimental data and simulated data.
    
    Parameters
    ----------
    data_arrays : array or list of 1D arrays
        Experimental data.
    sim_arrays : array or list of 1D arrays
        Simulated data. Must have same shape as `data_arrays`. Comparison will be
        made item by item.
    norm : boolean, optional
        Determines whether the comparison between data and simulation will be
        between normalized or un-normalized arrays. Default is True.
    absolute : boolean, optional
        Determines whether the function will return absolute differences (True)
        or the differences divided by the values of `data_arrays` (False). 
        Default is False.
    comparison : string
        Must be 'linear' or 'log'. Determines whether the function returns the
        difference in logs or the absolute difference between the two sets
        of arrays.
        
    Returns
    -------
    diffs : array or list of 1D arrays
        Returns the element-by-element relative differences between `data_arrays`
        and `sim_arrays` in the form of an array of the same shape as the input
        arrays. 
    """
    if np.all(sim_arrays==0):
        diffs = np.ones(sim_arrays.shape)*100000
        return diffs
    
    alert=False
    if norm:
        for i, array in enumerate(data_arrays):
            data_arrays = kin_kit.normalized(data_arrays, alert)
            sim_arrays = kin_kit.normalized(sim_arrays, alert)

    if absolute==False:
        diffs = abs( data_arrays[np.where(data_arrays !=0)] 
        - sim_arrays[np.where(data_arrays !=0)] ) / data_arrays[np.where(data_arrays !=0)]
        
    else:
        if comparison == 'linear':
            diffs = abs( data_arrays[np.where(data_arrays !=0)] - sim_arrays[np.where(data_arrays !=0)] ) 
        elif comparison == 'log':
            diffs = abs(np.log10(data_arrays[np.where(data_arrays > 0)]) 
            - np.log10(sim_arrays[np.where(data_arrays > 0)]))
        else: 
            print('Comparison must be linear or log.')
            return
    return diffs

def slice_by_time(arrays, timearray, lowlim=None, hilim=None):
    """
    Slices a set of arrays along values of a reference array with the same
    'long' dimension.
    """
    
    sliced_arrays = arrays.copy()
    
    lowidx, lowval = kin_kit.find_nearest(timearray, lowlim)
    hiidx, hival = kin_kit.find_nearest(timearray, hilim)
    
    if lowidx > hiidx:
        sliced_arrays = np.roll(sliced_arrays, -lowidx, axis=-1)
        hiidx += len(timearray) - lowidx
        if len(sliced_arrays.shape) != 1:
            return sliced_arrays[:, 0:hiidx+1]
        else:
            return sliced_arrays[0:hiidx+1]
    
    else:
        if len(sliced_arrays.shape) != 1:
            return sliced_arrays[:, lowidx:hiidx+1]
        else:
            return sliced_arrays[lowidx:hiidx+1]


def simulate_and_compare(varparams, varparamkeys, system, data_arrays, to, 
                         light, powers=None, irf_fwhm=40*ps, N_coarse=500, 
                         roll_value=0, comparison='linear', absolute=True, 
                         limits = None, norm=True, roll_criterion='max', 
                         maxavgnum=10, condensed_output=True):
    """
    Returns an array or list of arrays of differences between a set of
    simulated data and a set of experimental data, after aligning them. 
    The parameters and output of this functions are eventually to be passed to
    an optimization function such as ``scipy.optimize.least_sqares``.
    
    Parameters
    ----------
    varparams : list
        Rate equation parameters or other simulation parameters that determine
        the output of `system`. It does not need to be the comprehensive list 
        of parameters that define the system, simply the parameters with 
        respect to which the system output is to be optimized.
        varparams must have been obtained as list(dictionary.keys()) of the
        same dictionary that varparamkeys are extracted.
    varparamkeys : dictionary keys
        Names of varied parameters.
    system : object
        A system object to which `params` are passed, which contains a function
        that can output an array of PL signals. The system object needs to have
        been created previously in the code with dummy or guess parameters. It
        is then redefined with different parameters every time the 
        ``simulate_and_compare`` function is called.
    data_list : array or list of 1D arrays
        Array or list of 1D arrays containing experimental data to be compared 
        to the simulation. The output of ``system.PLsig`` must have the same 
        shape as `data_list`.
    to : dictionary
        Dictionary with time parameters. 
    light : object
        Excitation object that determines simulation
    powers : dictionary
        List of pulse powers at which experiment is conducted
    irf_fwhm : float, optional
        Width of Instrument Response Function used for convolution. See
        ``convolve_irf`` function. Default is 40 ps.
    N_coarse : integer
        N_coarse parameter of refined_simulation function.
    roll_value : float, optional
        Value on ``to['array']`` to which maxima of aligned arrays will be 
        shifted to.
    absolute : boolean, optional
        see ``elementwise_diff`` function. Default is False.
    limits : list of 2 values, optional
        lower and upper limit on time axis for which fitting is desired.
        Default is None (whole time axis fitted).
    norm : Boolean
        Determines whether each data/simulation transient should be normalized 
        with respect to itself. If false, all simulation and data traces are
        normalized by the simulation and data trace corresponding 
        to the highest power, respectively.
    roll_criterion : string, 'steep' or 'max'
        Determines whether the alignment of the two arrays happens by their 
        maximum or their steepest point. Default is 'max'
    maxavgnum : integer
        If `roll_criterion` is `"max"`, maxavgnum determines the `avgnum` 
        keyword in ``KinetiKit.kit.align_by_max().``
        
    Returns
    -------
    diffs : array or list of 1D arrays
        Array(s) of the differences or relative differences between the passed
        experimental data and simulation array(s).
    al_data_list : array or list of 1D arrays
        Aligned data arrays.
    al_sim_list : array or list of 1D arrays
        Aligned simulation arrays.
    """
      
    param_dict = kin_kit.dict_from_list(varparams, varparamkeys)
    system.update(**param_dict)
    if light is not None:
        pulse = light.pulse
    dtime = to['array'][::to['subsample']]
    
    global counter
    
    if counter%settings['display_counter_every'] == 0:
        if settings['display_counter']:
            print(counter)
    counter+=1
    
    if system.populations is None:
        pl, converged = sim.lib.simulate_func(system, dtime)
        sim_arrays = sim.lib.convolve_irf(pl, dtime, 
                                          fwhm=irf_fwhm)
    
    else:
        if powers is None:
    
            transient, converged = sim.lib.refined_simulation(system, to, light,
                                                          N_coarse=N_coarse)
            pl = system.PLsig(transient)
            sim_arrays = sim.lib.convolve_irf(pl, dtime, 
                                          fwhm=irf_fwhm)
            #sim_arrays /= max(sim_arrays)
            #data_arrays /= max(data_arrays)
            
        else:
            for i, power in enumerate(powers):
                pulse['power'] = power
                light = sim.lib.Excitation(pulse=pulse)
                transient, converged = sim.lib.refined_simulation(system, to, light,
                                                          N_coarse=N_coarse)
                pl_at_this_power = system.PLsig(transient)
                
                if i == 0:
                    pl = pl_at_this_power
                else:
                    pl = np.vstack((pl, pl_at_this_power))
                
            sim_arrays = sim.lib.convolve_irf(pl, dtime, 
                                          fwhm=irf_fwhm)    
    
    
    if roll_criterion == 'max':
        al_data_arrays = kin_kit.align_by_max(data_arrays, 
                                           dtime, 
                                           avgnum = maxavgnum,
                                           value = roll_value)
        
        al_sim_arrays = kin_kit.align_by_max(sim_arrays, 
                                           dtime,
                                           avgnum = maxavgnum,
                                           value = roll_value)
    
    elif roll_criterion == 'steep':
        al_data_arrays = kin_kit.align_by_steep(data_arrays,
                                                 dtime,
                                                 avgnum = maxavgnum,
                                                 value=roll_value)
        
        al_sim_arrays = kin_kit.align_by_steep(sim_arrays,
                                                dtime,
                                                avgnum = maxavgnum,
                                                value=roll_value)
        
    else: 
        print('Roll_criterion must be max or steep.')
    
    sim_arrays = kin_kit.make_2d(sim_arrays)
    data_arrays = kin_kit.make_2d(data_arrays)
    
    # divide all traces by maximum value of highest-power trace 
    sim_arrays /= max(sim_arrays[np.argmax(powers)]) 
    data_arrays /= max(data_arrays[np.argmax(powers)])  
    
    if limits is not None:    
        al_data_arrays = slice_by_time(al_data_arrays, dtime, limits[0], limits[1])
        al_sim_arrays = slice_by_time(al_sim_arrays, dtime, limits[0], limits[1])
    
    diffs = elementwise_diff(al_data_arrays, al_sim_arrays, 
                                 norm = norm, comparison=comparison,
                                 absolute = absolute)
        
    #print("diff : %0.3e"%(np.sum(diffs**2)/(to['N']/to['subsample'])))
    if condensed_output:
        return np.sum(diffs**2)
    else:
        return diffs

def sac_args(varparamkeys, system, data_arrays, to, 
                         light, powers=None, irf_fwhm=40*ps, N_coarse=500, roll_value=0, 
                         comparison='linear', absolute=True, limits = None,
                         norm=True, roll_criterion='max', maxavgnum=10,
                         condensed_output=True):
    """
    Returns a list of all but the first argument of ``simulate_and_compare()``.
    Used inside a minimization function like ``differential_evolution``, which
    does not accept keyword arguments, to selectively change some of the 
    arguments.
    """
    
    return varparamkeys, system, data_arrays, to, light, powers,  irf_fwhm, \
N_coarse, roll_value, comparison, absolute, limits, norm, roll_criterion, \
maxavgnum, condensed_output

def fit_leastsq(function, p0, args):
    # original idea by https://stackoverflow.com/a/21844726
    errfunc = function(p0, *args)
    
    pfit, pcov, infodict, errmsg, success = \
        sp.optimize.leastsq(function, p0, args=args, \
                          full_output=1)

    if (len(errfunc) > len(p0)) and pcov is not None:
        s_sq = (function(pfit, *args)**2).sum()/(len(errfunc)-len(p0))
        pcov = pcov * s_sq
    else:
        pcov = np.inf

    error = [] 
    for i in range(len(pfit)):
        try:
          error.append(np.absolute(pcov[i][i])**0.5)
        except:
          error.append( 0.00 )
    pfit_leastsq = pfit
    perr_leastsq = np.array(error) 
    return pfit_leastsq, perr_leastsq, infodict