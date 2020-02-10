
__all__ = ['simulate_func']

def simulate_func(function, t):
    """
    Simulates a system based on a mathematical function along a time axis.
    
    Required Parameters
    ----------
    function : system instance
        requires `PLsig` method
    t_obj : dict
        Key value combinations from sim.time module
    
    """    
    out = function.PLsig(t)
    
    return out, True
    