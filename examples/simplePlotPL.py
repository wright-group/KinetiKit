# -*- coding: utf-8 -*-
"""
Created on Fri Sep 21 15:12:50 2018

@author: Natalia Spitha
Displays a steady-state PL spectrum or a series thereof.
"""
import os

import numpy as np
import matplotlib 
import matplotlib.pyplot as plt

from KinetiKit.artists import plotparams

matplotlib.rcParams['agg.path.chunksize'] = 10000

dir_path = os.path.dirname(os.path.realpath(__file__))
#dir_path = r'C:\Users\...'
subfolder = r'ex_data'

plt.figure()
filenames = [r'1D_PL_a.asc', '1D_PL_b.asc']
plt.title(filenames[0].split('_PL')[0])
labels = ['n=1', 'n=2']
colors = ['blue','green']

for i, fname in enumerate(filenames):
    p = os.path.join(dir_path, subfolder, fname)
    data = np.genfromtxt(p, delimiter = ',', skip_footer = 22, unpack = True)
    wavelengths = data[0]
    pl = data[1]
    pl -= min(pl)
    pl /= max(pl)
    energies = 1240/wavelengths
    plt.plot(wavelengths,pl,label=labels[i],color=colors[i],alpha=0.5)
plt.xlim(480,780)
plt.xlabel('Wavelength (nm)')
plt.ylabel('PL intensity (a.u.)')
plt.legend(loc=1)
plt.tight_layout()
plt.show()