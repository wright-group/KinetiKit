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
        to=sim.time.linear(), N_coarse=500, power = 1e-6, irf_whm = 50*ps,
        data=None, power_list=[1], power_unit='microWatt', 
        align_by = 'steep', avgnum= 5, slidepower=False):
    
    args = p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14,\
        p15, p16, p17, p18, p19, p20,
    param_names = system.params().keys()
    trunc_args = args[:len(list(param_names))]
    params = kin_kit.dict_from_list(trunc_args, param_names)
    system.update(**params) # system is updated according to the parameters provided as *args
    
    if isinstance(power_list, np.float) or isinstance(power_list, np.int):
        power_list = [power_list]
        
    #--- Creating Time Array
    dtime = to['array'][::to['subsample']]
    
    if slidepower:
        power *= units[power_unit]
        light = sim.lib.Excitation(pulse={'power' : power})
        transient, converged = sim.lib.refined_simulation(system, to, light,
                                                  N_coarse=N_coarse)
        pl = system.PLsig(transient)
    
    else:
        for i, power in enumerate(power_list):
            power *= units[power_unit]
            light = sim.lib.Excitation(pulse={'power' : power})
            transient, converged = sim.lib.refined_simulation(system, to, light,
                                                      N_coarse=N_coarse)
            pl_at_this_power = system.PLsig(transient)
            
            if i == 0:
                pl = pl_at_this_power
            else:
                pl = np.vstack((pl, pl_at_this_power))
                
    sims = sim.lib.convolve_irf(pl, dtime, fwhm=irf_fwhm)   
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
                         mlabel='data_%0.3f'%(power_list[i])
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
                                 mlabel='sim_%0.3f'%(power_list[i])
                                 )
            except IndexError:
                sims_ended = True
                pass
            
            if sims_ended and data_ended:
                break
            
    return

def HeteroViz(system, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14,
        p15, p16, p17, p18, p19, p20,
        to=sim.time.linear(), N_coarse=500, power = 1e-6, irf_whm = 50*ps, 
        data=None, power_unit='microWatt', ids = ['layer 1', 'layer 2'],
        align_by = 'steep', avgnum= 5):
    
    args = p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14,\
        p15, p16, p17, p18, p19, p20,
    param_names = system.params().keys()
    trunc_args = args[:len(list(param_names))]
    params = kin_kit.dict_from_list(trunc_args, param_names)
    system.update(**params) # system is updated according to the parameters provided as *args
            
    #--- Creating Time Array
    dtime = to['array'][::to['subsample']]
    
    power *= units[power_unit]
    
    light = sim.lib.Excitation(pulse={'power' : power})
    transient, converged = sim.lib.refined_simulation(system, to, light,
                                                  N_coarse=N_coarse)
    pl = kin_kit.make_2d(system.PLsig(transient))
    sims = sim.lib.convolve_irf(pl, dtime, fwhm=irf_fwhm)   
    
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
                         mlabel='data_%s'%(ids[i])
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
                                 mlabel='sim_%s'%(ids[i])
                                 )
            except IndexError:
                sims_ended = True
                pass
            
            if sims_ended and data_ended:
                break
            
    return

def FuncViz(system, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14,
        p15, p16, p17, p18, p19, p20,
        to=sim.time.linear(), N_coarse=500, power = 1e-6, irf_fwhm=50*ps,
        data=None, power_list=[1], power_unit='microWatt', 
        align_by = 'steep', avgnum= 5, slidepower=False):
    
    args = p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14,\
        p15, p16, p17, p18, p19, p20,
    param_names = system.params().keys()
    trunc_args = args[:len(list(param_names))]
    params = kin_kit.dict_from_list(trunc_args, param_names)
    system.update(**params) # system is updated according to the parameters provided as *args
    
            
    #--- Creating Time Array
    dtime = to['array'][::to['subsample']]
    
    pl = system.PLsig(dtime)
    
    sims = sim.lib.convolve_irf(pl, dtime, fwhm=irf_fwhm)  
    
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
                         mlabel='data_%0.3f'%(power_list[i])
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
                                 mlabel='sim_%0.3f'%(power_list[i])
                                 )
            except IndexError:
                sims_ended = True
                pass
            
            if sims_ended and data_ended:
                break
            
    return
