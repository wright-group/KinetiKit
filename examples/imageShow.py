# -*- coding: utf-8 -*-
"""
Created on Tue Oct 17 13:11:32 2017

@author: Natalia Spitha

Script for displaying monochromatic .asc images from ANDOR solis in the 
correct format, with a spatial scale. Can be updated with the length/px factors
for new, calibrated imaging systems

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sps
from KinetiKit.artists import plotparams

dir_path = os.path.dirname(os.path.realpath(__file__))
subfolder = r'ex_data'
filename = 'monochr_image.asc'
spectral = False # is this a spatially-resolved PL spectrum or an optical image?

file = os.path.join(dir_path, subfolder, filename)

title = 'sample image'

magnifications = { 
#   'magn': (<vertical μm/px>, <horizontal μm/px>)
    '100x': (0.06027, 0.06035),
    '20x' : (0.29985,0.33418),
                 }

vcf, hcf = magnifications['20x']

data = np.genfromtxt(file, delimiter=',', skip_header=0, skip_footer=33, 
                     unpack=True)[:-1]

x = data[0]
vpixels = np.arange(0, data.shape[0]-1, 1)
vdists = vpixels * vcf #determined vertical μm/px
vd_start = vdists[0]; vd_end = vdists[-1]
signals = data[1:]

if spectral:
    threshold = np.max(signals)
    minimum= 100
    wavelengths = x
    wl_start = wavelengths[0]
    wl_end = wavelengths[-1]
    extent1 = [wl_start, wl_end, -vd_end/2, vd_end/2]
    x_label = r"Wavelength ($nm$)"
    signals[signals > threshold] = threshold
    signals[signals < minimum] = minimum
    signals -= minimum
    plt.xticks

else:
    threshold= np.max(signals)
    minimum = 100
    hpixels = x-1
    hpx_start = hpixels[0]
    hpex_end = hpixels[-1]
    hdists = hpixels * hcf #determined horizontal μm/px
    hd_start = hdists[0]; hd_end = hdists[-1]
    extent1 = [-hd_end/2,hd_end/2,-vd_end/2,vd_end/2]
    x_label = r"Distance ($\mu m$)"
    signals[signals > threshold] = threshold
    signals[signals < minimum] = minimum
    signals -= minimum
    

plt.figure()
plt.title(title)
plt.imshow(signals/np.max(signals),aspect='auto',cmap='bone', extent=extent1)
plt.colorbar()
plt.grid() #used to turn OFF the grid
if spectral:
    plt.grid(axis = 'x',color='lightskyblue',alpha=0.6,linestyle = '--')
plt.xlabel(x_label)
plt.ylabel(r"Distance ($\mu m$)")
plt.tight_layout()
plt.show()
