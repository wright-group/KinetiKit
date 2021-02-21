import numpy as np
from scipy.signal import convolve, savgol_filter
from matplotlib import pyplot as plt

from KinetiKit.units import units, MHz, fs, nm, uW, ps
import KinetiKit.kit as kin_kit

__all__ = ['convolve_irf']

# --- IRF Functions
def convolve_irf(pl, t, args):
    """
    Convolves simulated function with a specified Instrument Response Function.
    IRF can be either simulated as a simple Gaussian or as an exponentially 
    modified Gaussian with a diffusion tail. See ``build_irf`` for the default 
    values for IRF construction.
    
    Parameters
    ----------
    pl : numpy array
        The signal to be convolved must be a 1-D or 2-D array, with the lowest
        dimension having the same length as the time array`t`
    t : numpy array
        1D array representing the time axis.
    args : dictionary
        Arguments to be input into `build_irf` function. See definition of
        ``build_irf`` for default values.
    
    Returns
    -------
    An array of the shape of `pl` representing the signal convolved with the 
    IRF.
    """
    
    irf = build_irf(t, **args)
    if pl.ndim == 2:
        irf = irf[np.newaxis]
    return convolve_cyclic_boundary(pl, irf)


def build_irf(t, irf_type='Gauss', weighted=True, fwhm=55 * ps, tau_wt=40 * ps, 
              tau=650 * ps, b=0.13, c=0):
    """
    Constructs a Gaussian-like curve simulating an instrument response function
    of a time-resolved measurement system. 
    
    Parameters
    ----------
    t : numpy array
        The time array on which to build the IRF
    irf_type : string
        Type of IRF to be constructed; ```irf_type = 'Gauss'``` will return a 
        normal Gaussian curve, while ```irf_type = 'GaussDiff'``` will return
        an exponentially modified Gaussian with a diffusion tail.
    weighted : Boolean
        Whether to use a pure Gaussian or exponentially weighted Gaussian for 
        the main "peak" of the IRF. Default is True. 
        See en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution. 
    fwhm : float
        The FWHM that determines the variance of either type of curve. Must 
        have the same time units as `t`.
    tau_wt : float
        lifetime of decay of exponentially weighted Gaussian. Ignored if 
        `weighted` is set to False.
    tau : float
        lifetime of diffusion tail in the time units of `t`. Ignored if 
        `irf_type` is set to `'Gauss'`.
    b : float
        Relative contribution of diffusion tail to IRF. Ignored if 
        `irf_type` is set to `'Gauss'`.
    c : float
        Offset for IRF, recommended as zero.
    
        
    
    
    """
    
    tirf = t[t <= 3 * fwhm]
    tirf -= tirf.mean()
    if irf_type == 'Gauss':
        irf = kin_kit.Gauss(tirf, 1, 0, fwhm)
    elif irf_type == 'GaussDiff':
        irf = kin_kit.GaussDiff(t, np.max(t)/2, 1, fwhm=fwhm, tau=tau, 
                                b=b, c=c, weighted=weighted, tau_wt=tau_wt)
        window_length = int(0.006*len(t))
        if window_length%2 == 0:
            window_length += 1
        irf = savgol_filter(irf, window_length=window_length, polyorder=1)
    return irf / irf.sum()

def convolve_cyclic_boundary(data, irf):
    m = irf.size
    out = np.concatenate(
        (data[..., -m:], data[..., :], data[..., :m]), axis=-1)
    out = convolve(out, irf, mode='same')
    return out[..., m:-m]
