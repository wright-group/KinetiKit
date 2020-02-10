"""
Interconversion between arrays, lists, dictionaries
"""

import numpy as np

def dict_from_list(arglist, dictionary):
    """
    Returns a new dictionary with the same keys as `dictionary` where values
    are set from arglist. 
    
    Assumes that `arglist` has the same length as `dictionary`, and that the
    list items are ordered in the same order that `dictionary.keys()` were 
    defined. Alternatively, `dictionary` can be a dictionary.keys() instance.
    
    Parameters
    ----------
    arglist : array-like
        Must have the same number of elements as `dictionary`.
    dictionary : dictionary, or dictionary.keys()
        The output of this function will be a new dictionary with the same 
        set of keys as `dictionary`.
        
    Raises
    ------
    TypeEerror : If `dictionary` is not of the correct type.
    
    Returns
    -------
    new_dict : dictionary
        A new dictionary with same set of keys as `dictionary` and corresponding
        values set to values of `arglist`.
    """
    import collections
    new_dict = {}
    
    if isinstance(dictionary, dict):
        for i, key in enumerate(list(dictionary.keys())):
            new_dict[key] = arglist[i]
    elif isinstance(dictionary,collections.abc.KeysView):
        for i, key in enumerate(list(dictionary)):
            new_dict[key] = arglist[i]
    else:
        raise TypeError('dictionary must be of type dict or dict_keys.')
        #issue_error
    
    return new_dict

def list_to_array(listorarray):
    if isinstance(listorarray, list):
        return np.array(listorarray)
    elif isinstance(listorarray, np.ndarray):
        return listorarray
    #issue_error

 
def make_2d(array):
    array = list_to_array(array)
    if len(array.shape) ==1:
        return(np.reshape(array, (1,len(array))))
    elif len(array.shape) > 2:
        return np.concatenate(array)
    else:
        return array


