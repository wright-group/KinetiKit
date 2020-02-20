import numpy as np
from scipy.signal import convolve
from matplotlib import pyplot as plt

from KinetiKit.units import units, MHz, fs, nm, uW, ps, ns
import KinetiKit.kit as kin_kit

__all__ = ['convolve_irf']

# --- IRF Functions
def convolve_irf(pl, t, **kwargs):
    irf = build_irf(t, **kwargs)
    if pl.ndim == 2:
        irf = irf[np.newaxis]
    return convolve_cyclic_boundary(pl, irf)


def build_irf(t, fwhm=40 * ps,  diff=False, ampl=0.05, lifetime=1.2*ns):
    """
    

    Parameters
    ----------
    t : array
        Time array
    fwhm : float, optional
        The full width at half-maximum of the response function.
        The default is 40 * ps.
    diff : Boolean, optional
        Whether to include a diffusive tail in the IRF shape. The default is
        False.
    ampl : float, optional
        Number between 0 and 1. Amplitude of diffusive tail relative to peak
        of IRF. The default is 0. Ignored if ``diff=False``.
    lifetime : float, optional
        Lifetime, in seconds, of diffusive tail. The default is 0.85*ns. 
        Ignored if ``diff=False``.

    Returns
    -------
    IRF
        Shape of instrument response function, the area under which is 1.

    """
    sigma = fwhm/(2*np.sqrt(2*np.log(2)))
    tirf = t[t <= 3 * fwhm]
    tirf -= tirf.mean()
    irf = kin_kit.Gauss(tirf, 1, 0, fwhm)
    if diff:
        sigma = fwhm/(2*np.sqrt(2*np.log(2)))
        offset = 1.5*fwhm
        offset_idx = (np.abs(t-offset)).argmin()
        irf_gauss = kin_kit.Gauss(t, 1, offset, fwhm)
        lineshape = ampl*np.exp(-(t-offset)/lifetime)
        # quadratic solution for crossing point between lineshape and irf_gauss
        t_cross = 2*sigma**2/lifetime + 0.5*np.sqrt(4*sigma**4/lifetime**2 - 8*sigma**2*np.log(ampl))
        # plt.plot(t/ns, irf_gauss)
        # plt.plot(t/ns, lineshape)
        # print(t_cross/ns)
        irf_new = np.zeros(len(t))
        for i, irf_point in enumerate(irf_new):
            if t[i] <= t_cross+offset:
                irf_new[i] = irf_gauss[i]
            else:
                irf_new[i] = lineshape[i]
        irf = np.roll(irf_new, -offset_idx)
    return irf / irf.sum()


def convolve_cyclic_boundary(data, irf):
    m = irf.size
    out = np.concatenate(
        (data[..., -m:], data[..., :], data[..., :m]), axis=-1)
    out = convolve(out, irf, mode='same')
    return out[..., m:-m]

def plot_irf(time, fwhm=40 * ps,  diff=True, ampl=0.06, lifetime=0.85*ns, 
             norm=True, shift=1.6*ns):
    y = build_irf(time, fwhm=fwhm, diff=diff, ampl=ampl, lifetime=lifetime)
    if norm:
        y /= max(y)
    y = kin_kit.roll_by_array_shift(y, time, shift)
    plt.plot(time/ns, y, label = '%i ps, %0.2f, %.2f ns'%(fwhm/ps, ampl, lifetime/ns))
    plt.yscale('log')
    if norm:
        plt.ylim(1e-5, 2)
    plt.legend()    
