"""

Test for fitting a data object with a biexponential function. 
Requires numpy, scipy, matplotlib, and time packages along with TRPL package
"""

from KinetiKit import sim
from KinetiKit import artists as art
from KinetiKit import kit as kin_kit
from KinetiKit.units import ns, ps

#--- Creating Time Object
to = sim.time.linear(N=1000)
dtime = to['array'][::to['subsample']]

#--- Create system instance
system = sim.systems.Biexp()

#--- Parameters of simulation and initializing of system
params = {
           'A1': 1.9974e-1,
           'tau1': 1.774e-9,
           'tau2': 1e-10,
           'offset': 0,
                }

system.update(**params)

pl, converged = sim.lib.simulate_func(system, dtime)
sims = sim.lib.convolve_irf(pl, dtime, fwhm=40 * ps)

# Aligns simulation
aligned_sims = kin_kit.align_by_steep(sims, dtime, value = 0.5*ns)

#--- Plot
art.plot3scales.plot(dtime, aligned_sims, system, t_dict=to, ivtype='none', 
                     annotate=True, fig=None, ResetColorCyc=True, 
                     linewidth=6, opaq=0.4,
                     mlabel='biexp sim')
