"""
Displays a TRPL transient or a series thereof based on text files.
Author: Natalia Spitha, 6/3/2019
"""

import os
import numpy as np
import matplotlib.pyplot as plt

from KinetiKit.units import ns
from KinetiKit import data
from KinetiKit.artists import plotparams

norm = True

all_d = []
dir_path = os.path.dirname(os.path.realpath(__file__))
#dir_path = r'C:\Users...'
subfolder = r'ex_data'
filenames = ['pure_nometadata_TRPL.asc']
keys = ['perovskite']

for i, filename in enumerate(filenames):
    p = os.path.join(dir_path, subfolder, filename)
    all_d.append(data.lib.data_from_SPCM(p, key = keys[i], 
                                             skip_h=154, skip_f=884, coll=60))
plt.figure(figsize=(10,7))
ax = plt.subplot(111)
plt.title(keys[0])
for d in all_d:
    key = d.key
    time = d.x / ns
    y = d.y
    if norm:
        y /= max(y)
    plt.plot(time,y,label=key)
ax.set_xlim(0,12.5)
ax.set_xlabel("Time (ns)")
ax.set_ylabel("Counts")
ax.set_yscale('log')
plt.legend()
plt.show()


