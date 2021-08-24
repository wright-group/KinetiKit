import time

import numpy as np

from KinetiKit import sim


__all__ = ['simulate_until_steady', 'simulate', 
'simulate_for_cycles', 'simulate_until_steady',
'refined_simulation']


def simulate(p0, RE_set, t_obj, light_obj):
    """
    Simulates ONE pulse cycle given an initial population, rate equation, 
    time, and excitation parameters.
    """
    
    t = t_obj['array']
    dt = t_obj['dt']    
    subsample = t_obj['subsample']

    # --- Incident Photons --- #
    pulse = light_obj.gen_pulse(t[t <= light_obj.pulse_window], center=0.5*light_obj.pulse_window, stepsize=dt)
    cw = light_obj.cw_power

    # --- Population arrays --- #
    out = np.zeros((RE_set.popnum, t.size // subsample))
    p_current = p0.copy() if p0 is not None else np.zeros(RE_set.shape)
    
    for i in range(t.size):
        p_previous = p_current.copy()

        photons = cw
        if i < pulse.size:
            photons += pulse[i]
        if i<5:
            pass
            #print(i, p_current, RE_set.rate(p_previous, photons))
        p_current += dt * RE_set.rate(p_previous, photons)
        
        if i % subsample == 0:
            out[:, i // subsample] = p_current

    return out


def simulate_for_cycles(RE_set, t_obj, light_obj, numcycles=1, p0=None):
    """
    Performs a simulation given pulse and rate equation parameters up to a given
    number of pulse cycles. Notifies when steady state is reached but does not
    stop simulation until numcycles is reached. Designed for testing new simulation
    models for convergence.

    Required Parameters
    ----------
    RE_set : system instance
        requires `rate` method
    t_obj : dict
        Key value combinations from sim.time module
    light_obj : excitation object
        Contains information about the pulsed and CW excitation for the 
        simulated experiment.
    
    Optional Parameters
    ----------
    numcycles : integer
        Number of cycles for which simulation occurs. Default is 1.
    p0 : array 
        Initial population from which to start the simulation. Setting p0 to 
        None defaults to a population of zeros. p0 must have the length of
        `RE_set.popnum`.

    Returns
    ----------
    current : array (2D)
        The final iteration of simulation.  Shape of (npop, ntime).
    converged : bool
        Whether or not the simulation converged within numcycles loops

    """
    converged = False
    current = None
    if p0 is None:
        p0 = np.zeros(RE_set.popnum) 
    else:
        p0 = p0
    notified=False
    for c in range(1, numcycles+1):
        current = simulate(p0, RE_set, t_obj, light_obj)
        converged = isSteady(p0, current[:, -1])
        if converged and not notified:
            print("Reached steady after %i cycles. Continuing simulation"%c)
            notified = True
        p0 = current[:, -1]
    return current, converged

def simulate_until_steady(RE_set, t_obj, light_obj, p0=None, verbose=False):
    """
    Performs a simulation given pulse and rate equation parameters until steady
    state is reached.

    Required Parameters
    ----------
    RE_set : system instance
        requires `rate` method
    t_obj : dict
        Key value combinations from sim.time module
    light_obj : excitation object
        Contains information about the pulsed and CW excitation for the 
        simulated experiment.

    Optional Parameters
    ----------
    p0 : array-like
        initial carrier populations; must have the same length as 
        RE_set.populations. Default is all zeros.
    verbose : boolean
        Whether to print debug-friendly informative messages. Default is False.

    Returns
    ----------
    current : array (2D)
        The final iteration of simulation.  Shape of (npop, ntime).
    converged : bool
        Whether or not the simulation converged within light_obj.numCycles loops

    """
    converged = False
    toofast=False
    current = None
    if p0 is None:
        p0 = np.zeros(RE_set.popnum) 
    else:
        p0 = p0
    for c in range(1, light_obj.numcycles+1):
        current = simulate(p0, RE_set, t_obj, light_obj)
        if (current<0).any():
            #print('Too fast')
            toofast=True
            current = np.zeros(current.shape)-1
            break
        converged = isSteady(p0, current[:, -1], tol=0.001)
        if converged:
            if verbose:
                print("Reached steady after %i cycles"%c)
            break
        p0 = current[:, -1]
    if not converged and not toofast:
        print('Failed to reach steady state after %i cycles'%light_obj.numcycles)
    return current, converged

def refined_simulation(RE_set, t_obj, light_obj, N_coarse=500, doubleSearch=False, verbose=False):
    """
    Performs a simulation with coarse time step until steady state is reached, 
    and then performs a final simulation with a finer time step.

    Required Parameters
    ----------
    RE_set : system instance
        requires `rate` method
    t_obj : dict
        Key value combinations from sim.time module
    light_obj : excitation object
        Contains information about the pulsed and CW excitation for the 
        simulated experiment.
    N_coarse : integer
        Number of time steps per cycle taken during the coarse simulation.
        Default is 500.
        
    Optional Parameters
    -------------------
    doubleSearch : boolean
        Specifies whether the steady state population found by the coarse 
        simulation should be taken as is and only be integrated on once 
        (False), or whether a new search for steady state should be performed
        in the fine time scale (True). Setting `doubleSearch = True` guarantees
        that you will reach the same steady state as you would if you ran the 
        entire simulation at the fine time scale, but is unnecessary if the
        prediction of the steady state population by the coarse simulation is
        accurate. Default is False.
    verbose : boolean
        Whether to print debug-friendly informative messages. Default is False.
    
    Returns
    ----------
    fine_sim : array (2D)
        The final iteration of simulation.  
        Shape of (RE_set.popnum, `N_fine/t_obj['subsample']`).
    converged : bool
        Whether or not the simulation converged within t_obj.numCycles loops
    """
    
    #N_fine = t_obj['N']
    to_coarse = sim.time.update_linear(t_obj, **{'N' : N_coarse})
    tc_start = time.process_time()
    coarse_sim, converged = simulate_until_steady(RE_set, to_coarse, light_obj, verbose=verbose)
    coarse_p0 = coarse_sim.transpose()[np.argmin(coarse_sim[0])]
    tc_end = time.process_time()
    
    if (coarse_p0==-1).any(): # checks if coarse simulation led to "toofast" conditions
        #print('Did not refine simulation')
        return np.zeros((len(coarse_p0), t_obj['N'])), converged
    to_fine = t_obj
    
    tf_start = time.process_time()
    if doubleSearch:
        fine_sim, converged = simulate_until_steady(RE_set, to_fine, light_obj, p0=coarse_p0, verbose=verbose)
    else: 
        fine_sim = simulate(coarse_p0, RE_set, to_fine, light_obj)
    tf_end = time.process_time()
    
    #print('Coarse time (s):', tc_end-tc_start)
    #print('Fine time (s):', tf_end-tf_start)
    
    return fine_sim, converged


# --- Functions assisting integration
def isSteady(previous, current, tol=0.001):
    ids = current != 0
    dif = np.abs(current[ids] - previous[ids]) / current[ids]
    return np.all(dif < tol)

