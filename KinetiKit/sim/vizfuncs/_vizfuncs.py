# -*- coding: utf-8 -*-
"""
Functions that output an overlay of simulation and data given a system, 
a set of y data, and parameters. 

@author: Natalia Spitha
"""

import numpy as np

from KinetiKit import sim
from KinetiKit import artists as art
from KinetiKit.artists.plotparams import colors as color_range
from KinetiKit import kit as kin_kit
from KinetiKit.units import units, ns, ps


def MonoViz(system, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14,
        p15, p16, p17, p18, p19, p20,
        to=sim.time.linear(), N_coarse=500, power = 1e-6, irf_args = {},
        data=None, power_unit='microWatt', 
        align_by = 'steep', avgnum= 5, xmin=0.1, xmax=None, 
        ymin=1e-3, ymax=1.2):
    
    args = p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14,\
        p15, p16, p17, p18, p19, p20,
    param_names = system.params().keys()
    trunc_args = args[:len(list(param_names))]
    params = kin_kit.dict_from_list(trunc_args, param_names)
    system.update(**params) # system is updated according to the parameters provided as *args
           
    #--- Creating Time Array
    dtime = to['array'][::to['subsample']]
    
    if isinstance(power, np.float) or isinstance(power, np.int):
        power_list = [power]
        power *= units[power_unit]
        light = sim.lib.Excitation(pulse={'power' : power})
        transient, converged = sim.lib.refined_simulation(system, to, light,
                                                  N_coarse=N_coarse)
        pl = system.PLsig(transient)
    
    elif isinstance(power, list):
        power_list = power
        for i, power in enumerate(power):
            power *= units[power_unit]
            light = sim.lib.Excitation(pulse={'power' : power})
            transient, converged = sim.lib.refined_simulation(system, to, light,
                                                      N_coarse=N_coarse)
            pl_at_this_power = system.PLsig(transient)
            
            if i == 0:
                pl = pl_at_this_power
            else:
                pl = np.vstack((pl, pl_at_this_power))
                
    sims = sim.lib.convolve_irf(pl, dtime, irf_args)   
    # Aligns data with sim either by max. or steep
    if align_by == 'steep':
        aligned_sims = kin_kit.align_by_steep(sims, dtime, value = 0.5*ns, avgnum=avgnum)
    
    elif align_by == 'max':
        aligned_sims = kin_kit.align_by_max(sims, dtime, value = 0.5*ns, avgnum=avgnum)
    else:
        print('\'align_by\' must be \'steep\' or \'max\'')
    aligned_sims = kin_kit.make_2d(aligned_sims)   
    
    fig=None
    data_ended = False; sims_ended = False
    for i in range(20):
        if data is None:
            pass
        else:
            data=kin_kit.make_2d(data)
            try:
                d = data[i]
                fig=art.plot3scales.plot(dtime, d, sys_obj=None, t_dict=None,
                         annotate=False, fig=fig, linewidth=1, 
                         color=color_range[i],
                         mlabel='data_%0.3f'%(power_list[i]),
                         xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
                         
            except IndexError:
                data_ended = True
                pass
                
            try:
                aligned_sim = aligned_sims[i]
                art.plot3scales.plot(dtime, aligned_sim, system, t_dict=to, ivtype='none', 
                                 annotate=True, fig=fig, ResetColorCyc=i==0, 
                                 color=color_range[i],
                                 linewidth=8, opaq=0.4,
                                 mlabel='sim_%0.3f'%(power_list[i]),
                                 xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
                                 
            except IndexError:
                sims_ended = True
                pass
            
            if sims_ended and data_ended:
                break
            
    return

def HeteroViz(system, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14,
        p15, p16, p17, p18, p19, p20,
        to=sim.time.linear(), N_coarse=500, power = 1e-6, irf_args = {}, 
        data=None, power_unit='microWatt', ids = ['layer 1', 'layer 2'],
        light=sim.lib.Excitation(),
        align_by = 'steep', avgnum= 5, xmin=0.1, xmax=None, ymin=1e-3, ymax=1.2):
    
    args = p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14,\
        p15, p16, p17, p18, p19, p20,
    param_names = system.params().keys()
    trunc_args = args[:len(list(param_names))]
    params = kin_kit.dict_from_list(trunc_args, param_names)
    system.update(**params) # system is updated according to the parameters provided as *args    
    #--- Creating Time Array
    dtime = to['array'][::to['subsample']]
    
    power *= units[power_unit]
    light = light.updated_with(pulse={'power' : power})
    transient, converged = sim.lib.refined_simulation(system, to, light,
                                                  N_coarse=N_coarse)
    pl = kin_kit.make_2d(system.PLsig(transient))
    sims = sim.lib.convolve_irf(pl, dtime, irf_args)   
    
    # Aligns data with sim either by max. or steep
    if align_by == 'steep':
        aligned_sims = kin_kit.align_by_steep(sims, dtime, value = 0.5*ns, avgnum=avgnum)
    
    elif align_by == 'max':
        aligned_sims = kin_kit.align_by_max(sims, dtime, value = 0.5*ns, avgnum=avgnum)
    else:
        print('\'align_by\' must be \'steep\' or \'max\'')
    
    aligned_sims = kin_kit.make_2d(aligned_sims)   
    
    fig=None
    
    data_ended = False; sims_ended = False
    for i in range(20):
        if data is None:
            pass
        else:
            data=kin_kit.make_2d(data)
            try:
                d = data[i]
                fig=art.plot3scales.plot(dtime, d, sys_obj=None, t_dict=None,
                         annotate=False, fig=fig, linewidth=1, 
                         color=color_range[i],
                         mlabel='data_%s'%(ids[i]),
                         xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
            except IndexError:
                data_ended = True
                pass
                
            try:
                aligned_sim = aligned_sims[i]
                art.plot3scales.plot(dtime, aligned_sim, system, t_dict=to, ivtype='none', 
                                 annotate=True, fig=fig, ResetColorCyc=i==0, 
                                 color=color_range[i],
                                 linewidth=8, opaq=0.4, 
                                 mlabel='sim_%s'%(ids[i]),
                                 xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
                
            except IndexError:
                sims_ended = True
                pass
            
            if sims_ended and data_ended:
                break
            
    return

def MultiPowerViz(system, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14,
        p15, p16, p17, p18, p19, p20,
        to=sim.time.linear(), N_coarse=500, pulse_power = 1e-6, cw_power=0, irf_args = {'fwhm':55*ps}, 
        data=None, power_unit='microWatt', ids = ['layer 1', 'layer 2'], 
        light=sim.lib.Excitation(),
        align_by = 'steep', avgnum= 5, xmin=0.1, xmax=None, ymin=1e-3, ymax=1.2):
    """
    more functionality to allow for CW and pulsed power variations
    """
    args = p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14,\
        p15, p16, p17, p18, p19, p20,
    param_names = system.params().keys()
    trunc_args = args[:len(list(param_names))]
    params = kin_kit.dict_from_list(trunc_args, param_names)
    system.update(**params) # system is updated according to the parameters provided as *args
            
    #--- Creating Time Array
    dtime = to['array'][::to['subsample']]
    
    pulse_power *= units[power_unit]
    cw_power *= units[power_unit]
    light = light.updated_with(pulse={'power' : pulse_power}, cw = {'power': cw_power})
    transient, converged = sim.lib.refined_simulation(system, to, light,
                                                  N_coarse=N_coarse)
    pl = kin_kit.make_2d(system.PLsig(transient))
    sims = sim.lib.convolve_irf(pl, dtime, args=irf_args)   
    
    # Aligns data with sim either by max. or steep
    if align_by == 'steep':
        aligned_sims = kin_kit.align_by_steep(sims, dtime, value = 0.5*ns, avgnum=avgnum)
    
    elif align_by == 'max':
        aligned_sims = kin_kit.align_by_max(sims, dtime, value = 0.5*ns, avgnum=avgnum)
    else:
        print('\'align_by\' must be \'steep\' or \'max\'')
    
    aligned_sims = kin_kit.make_2d(aligned_sims)   
    
    fig=None
    
    data_ended = False; sims_ended = False
    for i in range(20):
        if data is None:
            pass
        else:
            data=kin_kit.make_2d(data)
            try:
                d = data[i]
                fig=art.plot3scales.plot(dtime, d, sys_obj=None, t_dict=None,
                         annotate=False, fig=fig, linewidth=1, 
                         color=color_range[i],
                         mlabel='data_%s'%(ids[i]),
                         xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
            except IndexError:
                data_ended = True
                pass
                
            try:
                aligned_sim = aligned_sims[i]
                art.plot3scales.plot(dtime, aligned_sim, system, t_dict=to, ivtype='none', 
                                 annotate=True, fig=fig, ResetColorCyc=i==0, 
                                 color=color_range[i],
                                 linewidth=8, opaq=0.4, 
                                 mlabel='sim_%s'%(ids[i]),
                                 xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
                
            except IndexError:
                sims_ended = True
                pass
            
            if sims_ended and data_ended:
                break
            
    return
    
def FuncViz(system, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14,
        p15, p16, p17, p18, p19, p20,
        to=sim.time.linear(), N_coarse=500, power = 1e-6, irf_args={},
        data=None, power_list=[1], power_unit='microWatt', 
        align_by = 'steep', avgnum= 5, slidepower=False,
        xmin=0.1, xmax=None, ymin=1e-3, ymax=1.2):
    
    args = p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14,\
        p15, p16, p17, p18, p19, p20,
    param_names = system.params().keys()
    trunc_args = args[:len(list(param_names))]
    params = kin_kit.dict_from_list(trunc_args, param_names)
    system.update(**params) # system is updated according to the parameters provided as *args
    
            
    #--- Creating Time Array
    dtime = to['array'][::to['subsample']]
    
    pl = system.PLsig(dtime)
    
    sims = sim.lib.convolve_irf(pl, dtime, irf_args)  
    
    # Aligns data with sim either by max. or steep
    if align_by == 'steep':
        aligned_sims = kin_kit.align_by_steep(sims, dtime, value = 0.5*ns, avgnum=avgnum)
    
    elif align_by == 'max':
        aligned_sims = kin_kit.align_by_max(sims, dtime, value = 0.5*ns, avgnum=avgnum)
    else:
        print('\'align_by\' must be \'steep\' or \'max\'')
    aligned_sims = kin_kit.make_2d(aligned_sims)   
    
    fig=None
    data_ended = False; sims_ended = False
    for i in range(20):
        if data is None:
            pass
        else:
            data=kin_kit.make_2d(data)
            try:
                d = data[i]
                fig=art.plot3scales.plot(dtime, d, sys_obj=None, t_dict=None,
                         annotate=False, fig=fig, linewidth=1, 
                         color=color_range[i],
                         mlabel='data_%0.3f'%(power_list[i]),
                         xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax
                         )
            except IndexError:
                data_ended = True
                pass
                
            try:
                aligned_sim = aligned_sims[i]
                art.plot3scales.plot(dtime, aligned_sim, system, t_dict=to, ivtype='none', 
                                 annotate=True, fig=fig, ResetColorCyc=i==0, 
                                 color=color_range[i],
                                 linewidth=8, opaq=0.4,
                                 mlabel='sim_%0.3f'%(power_list[i]),
                                 xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax
                                 )
            except IndexError:
                sims_ended = True
                pass
            
            if sims_ended and data_ended:
                break
            
    return
