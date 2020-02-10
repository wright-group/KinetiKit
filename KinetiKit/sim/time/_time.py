from KinetiKit.units import ns, ps, fs, MHz
import numpy as np


def linear(period=12.5 * ns,N=5000,
    subsample=1, pulse_window=2 * ps):
    """Time parameters for simulation.
    Assumes linearly spaced datapoints.

    Optional Parameters
    ----------
    period : float
        Duration of the time axis, in units of seconds.  For pulsed
        simulations, repetition rate is `1 / period`.  Default value is
        12.5 ns.
    N : int
        Number of time points sampled across `period`. Default is 50000.
    subsample : int
        How many time points are simulated between recorded data points.
        Default is 20. Simulation output will have `N // subsample` time
        points.
    pulse_window : float
        Defines the interval of the time array over which the pulse must 
        be accounted for. Units of seconds. Default value is 40 ps.

    Returns
    -------
    dictionary
        A collection of parameters relating to numerical integration and
        simulation output.
    """
    dic = {
        'period': period,
        'N': N,
        'subsample': subsample
    }

    dic['dt'] = period / N
    dic['stepsize'] = dic['dt'] * subsample
    dic['array'] = np.linspace(0, period, N + 1)[:-1]

    return dic

def get_linear_args(timeobject):
    timekeys = ['period', 'N', 'subsample']
    linear_args = {}
    for key, val in timeobject.items():
        if key in timekeys:
            linear_args[key] = val
    return linear_args

def update_linear(timeobject, **kwargs):
    timekeys = ['period', 'N', 'subsample']
    new_args = get_linear_args(timeobject)
    for key, val in kwargs.items():
        if key in timekeys:
            new_args[key] = val
        else:
            print("Expected one of the following: \n", timekeys)
    return linear(**new_args)


"""
Snippet to calculate computation time
Copy and paste:
    
time_start = time.clock()

at start of process and:
    
time_elapsed = (time.clock() - time_start)

at end of program.
"""