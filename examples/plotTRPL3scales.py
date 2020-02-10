"""
Displays a TRPL transient or a series thereof on 3 scales.
One may play around with the plot colors to their satisfaction.
Author: Natalia Spitha, 6/3/2019
"""

import os

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import KinetiKit.kit as kin_kit
from KinetiKit.units import ns
from KinetiKit import data, sim
import KinetiKit.artists as art

matplotlib.rcParams['lines.linewidth'] = 1

#--- Define time dictionary
to = sim.time.linear(period=12.5*ns, N=1000)
dtime = to['array']

all_data = []
all_y = []
dir_path = os.path.dirname(os.path.realpath(__file__))
#dir_path = r'C:\Users\...'
subfolder = 'ex_data'
destfolder = os.path.join(dir_path, 'ex_output')
commonkey = 'plotTRPL3scales'

filenames = ['het_n=1_0.54uW_TRPL.asc',
             'het_n=2_0.54uW_TRPL.asc', 
             'het_n=1_1.00uW_TRPL.asc', 
             'het_n=2_1.00uW_TRPL.asc',
             ]

colors = ['#E24A33', '#348ABD', '#988ED5', '#8EBA42', '#FBC15E', '#777777', '#FFB5B8']

for filename in filenames:
    key = filename.split('_TRPL')[0] #can be anything you want
    p = os.path.join(dir_path, subfolder, 'hetero', filename)
    do = data.lib.data_from_SPCM(p, key=key, 
                                  skip_h=154, skip_f=884)
    do.interp(dtime)
    all_y.append(do.y)
    all_data.append(do)
    
all_y = np.array(all_y)
    
fig = None
for i, d in enumerate(all_data):
    y = all_y[i]; key = d.key
    y = kin_kit.roll_by_array_shift(y, dtime, 2.2*ns, direction='back')
    markers = [True, False][i%2]
    fig=art.plot3scales.plot(dtime, y, fig=fig, ivtype = 'none',
                         unit = 'ns', annotate=False, norm=True, 
                         mlabel=key, markers=markers, color=colors[i//2],
                         opaq=0.6)
    plt.tight_layout()

plt.savefig(destfolder+commonkey+'.png', format='png' )