import numpy as np
from scipy.signal import convolve

from KinetiKit.units import units, MHz, fs, nm, uW, ps
import KinetiKit.kit as kin_kit

__all__ = ['convolve_irf']

# --- IRF Functions
def convolve_irf(pl, t, fwhm=40 * ps):
    irf = build_irf(t, fwhm=fwhm)
    if pl.ndim == 2:
        irf = irf[np.newaxis]
    return convolve_cyclic_boundary(pl, irf)


def build_irf(t, fwhm=40 * ps):
    tirf = t[t <= 3 * fwhm]
    tirf -= tirf.mean()
    irf = kin_kit.Gauss(tirf, 1, 0, fwhm)
    return irf / irf.sum()


def convolve_cyclic_boundary(data, irf):
    m = irf.size
    out = np.concatenate(
        (data[..., -m:], data[..., :], data[..., :m]), axis=-1)
    out = convolve(out, irf, mode='same')
    return out[..., m:-m]
