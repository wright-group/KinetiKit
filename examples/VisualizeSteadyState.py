"""
Useful for visualizing how steady state is reached on a cycle-per-cycle basis.
Author: Natalia Spitha
"""
import time

import numpy as np
import matplotlib.pyplot as plt

from KinetiKit import sim, artists
from KinetiKit import kit as kin_kit

from KinetiKit.units import ns

plot_output = False

var_values = [1,2,3,4,5,6,7, 8, 9, 10, 11, 12]  # values that the varied parameter will take
fwhm = 100e-12  # irf width [seconds]

to = sim.time.linear(period=12.5*ns, N = 1000, subsample=1)
dtime = to['array'][::to['subsample']]

# --- Simulation Params
params = {
    'k_ann':2e8,
    'k_dis': 1e9,
    'k_rec': 1e4,
    'cs': 0.2,
}
system = sim.systems.MonoRecX(**params)
excitation = sim.lib.Excitation()


plt.figure()
ax = plt.subplot(111)
ax2 = ax.twinx()
alphas = np.linspace(0.1,1,num=len(var_values))
for i, var_val in enumerate(var_values):
    excitation = excitation.updated_with(numcycles=var_val)
    transient, converged = sim.lib.refined_simulation(system, to, 
                                                      excitation,
                                                      N_coarse=500,
                                                      doubleSearch=False)
    transient = kin_kit.roll_by_array_shift(transient, dtime, 0.3*ns)

    rpl = system.PLsig(transient)
    
    pl = sim.lib.convolve_irf(rpl, dtime, fwhm=fwhm)
    
    xs = transient[0]; es = transient[1]; hs = transient[2]
    
    xs /= 1e-14 #divide by excitation volume
    es /= 1e-14
    hs /= 1e-14

    alpha = alphas[i]
    ax.plot(dtime/ns, xs, color='blue', alpha=alpha, label=var_val)
    ax2.plot(dtime/ns, es, color='red', alpha = alpha,)

if plot_output:
    ax.plot(dtime/ns, kin_kit.normalized(rpl), ':', color='black', linewidth=2, alpha=0.7, label='PL_signal')
    ax.plot(dtime/ns, kin_kit.normalized(pl), color='orange', linewidth=8, alpha=0.4, label='convolved PL')

ax.legend(loc=1)
plt.title('Reaching Steady State Visualization')
ax.set_xlabel('Time (ns)')
ax.set_ylabel('Norm. X population', color='blue')
ax2.set_ylabel('Norm E/H population', color='red')
ax.set_yscale('log'); ax2.set_yscale('log'); ax2.set_xscale('log')
# ax.set_ylim(0.0, 1.1)
# ax2.set_ylim(0.0,1.1)
plt.xlim(0.2, 12.5)
plt.legend(title='cycle')
plt.tight_layout()
    
plt.show()
