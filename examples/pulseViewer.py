"""
Pulse Viewer

Views a pulse under a given time resolution to display its value in photons/second.
If N is increased sufficiently for several time steps to fit "inside" a pulse,
the pulse will have a Gaussian shape. Otherwise, the pulse is treated as a 
step function. In both cases, the area under the curve contains the same amount
of photons that a "real" pulse of the defined specifications would contain.

Author: Natalia Spitha
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from KinetiKit import sim
import KinetiKit.artists as art
from KinetiKit.units import ps, MHz, fs, nm

matplotlib.rcParams['axes.facecolor']='black'
matplotlib.rcParams['lines.linewidth']=4

light_obj = sim.lib.Excitation(pulse={'power':1e-6, 
                                      'fwhm' : 100*fs,
                                      'wavelength':400*nm, 
                                      })
t_dict = sim.time.linear(N=1000000, subsample=1)
t = t_dict['array']
dt = t_dict['dt']

if light_obj.pulse_window < dt:
    view_window = dt
else:
    view_window = light_obj.pulse_window

pulse_t = t[t <= view_window]

print('pulse_t', len(pulse_t))

pulse = light_obj.gen_pulse(t[t <= view_window], 
                            center=view_window/2, stepsize=dt )

print('total photons delivered = ',  np.sum(pulse)*dt)

plt.figure()
ax=plt.subplot(1,1,1)
plt.title('Pulse Viewer (%0.2f μW, %i fs, N= %i)'%(light_obj.pulse_energy*80*MHz*1e6, 
          light_obj.pulse_fwhm/fs, t_dict['N']))

if len(pulse)==1:
    plotted_time = np.array([-view_window, 0, view_window, 2*view_window])
    plotted_pulse = np.array([0, 0, pulse[0], 0])
    plt.step(plotted_time/ps, plotted_pulse, color='white', alpha=0.5)
    plt.ylim(-pulse/50, 2*pulse)
    plt.xlim(-0.2*view_window/ps, 2.1*view_window/ps)
    plt.fill_between(plotted_time/ps, np.zeros(4), plotted_pulse,
                     step='pre', color='aquamarine', alpha=0.8)
else: 
    plt.plot(pulse_t/ps, pulse, color='white', alpha = 0.5)
    plt.fill_between(pulse_t/ps, np.zeros(len(pulse)), pulse,
                     step='pre', color='aquamarine', alpha=0.8)
plt.text(0.8, 0.9, r'Photons Delivered:', ha='center', va='top', color='white',
         transform = ax.transAxes)
plt.text(0.8, 0.85, '%i'%(np.sum(pulse)*dt), ha='center', va='top', color='white',
         transform = ax.transAxes)
plt.xlabel('Time (ps)')
plt.ylabel('Incident Photons/s')
plt.tight_layout()
plt.show()