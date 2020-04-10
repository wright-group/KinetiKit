import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

from KinetiKit.sim.time import linear
from KinetiKit.artists import plotparams
from KinetiKit.units import units
import KinetiKit.kit as kin_kit


def plot(x, y, sys_obj=None, t_dict=None, other_params=None, ivtype='none', 
         ivkey='k_ann', unit='ns', offset=None, norm=True, 
         fig=None, annotate=True, ResetColorCyc=False, 
         opaq=1, IgnoreLabel=False, mlabel=None, color=None,
         markers=False, linewidth=None, xmin=None, xmax=None, ymin=None, ymax=None):
    """
    A versatile function for plotting series of data or simulations on three
    scale types: linear, log-linear, and log-log. It contains several arguments
    relating to simulation parameters that can help automatically annotate
    the output plot.
    
    Parameters
    ----------
    x : 1D array
        Time array of data or simulation
    y : 1D array
        Signal array of data or simulation. Must be of same size as x.
    sys_obj : sim.system object, optional
        Determines the set of parameters such as rate equation constants that 
        characterize the simulation. It is used for properly annotating the 
        graph if a simulation is displayed. Default is None.
    t_dict : dictionary, optional
        Determines the set of parameters related to time which affect the 
        simulation, e.g. `reprate`, `N`, and `subsample`. Default is
        `sim.time.linear()`. It is used to properly annotate the time 
        parameters of the simulation. In a future iteration of this function, 
        the argument `x` could become optional and be replacable by 
        `t_dict['array']` If None, then `t_dict` is later defined as 
        `sim.time.linear`.
    other_params : dictionary, optional
        If a parameter other than system or time parameters is to be varied 
        (e.g. emission color for a single system, laser power, etc.), this
        dictionary should contain the keys and values for the additional 
        parameters that are relevant.
    ivtype : a string out of `'eq'`, `'tm'`, `'other'`, and `'none'`, optional
        Specifies the type of parameter that is varied across different traces
        that are plotted. `'eq'` specifies a simulation parameter, `'tm'` a time 
        parameter, `'other'` an additional parameter, and `'none'` the lack of 
        any varied parameter. Default is 'none'. 
    ivkey : string representing an existing dictionary key, optional
        If `ivtype` is 'eq', 'tm', or 'other', ivkey specifies the key of the 
        parameter that is varied in the system, time, or other parameter 
        dictionary, respectively. Default is None. If any parameter is to be 
        varied, `ivtype` and `ivkey` must always be defined.
    unit : string representing key of ``TRPL.units.units`` dictionary, optional
        Determines the unit written on x-axis and the value by which to divide 
        time array. Default is 'ns'.
    offset : array or length of length 3, optional
        If desired, plots can be made with a time offset specified in
        appropriate time units. The values of `offset[0]`, `offset[1]`, and
        `offset[2]` determine the shifts in the linear, log-linear, and log-log
        plots, respectively. Default is None.
    norm : Boolean, optional
        Specifies whether normalized or unnormalized graphs should be plotted. 
        Default is True.
    fig : Figure object or None
        Determines whether the plot traces should be made on a new figure or on
        an existing set of axes. Default is None, in which case a new figure 
        is plotted. If a figure object is provided, it must have (1, 3) 
        subplots
    annotate : Boolean, optional
        Specifies whether to display the simulation and time parameters at the
        bottom of the three graphs. Default is True. The size of the figure 
        depends on the existence or absence of annotation.
    ResetColorCyc : Boolean, optional
        If True, the color cycler will be reset to the first value specified 
        in ``artists.plotparams``. Should be used in conjunction with varying
        opaq so that traces with similar color can always be distinguished. 
        Can be useful when plotting two traces (e.g. emission color of the 
        same system) and want to show how each of them varies with a given 
        parameters. See ``simRates_hetero_vary.py`` for an example.
    opaq : Float between 0 and 1, optional
        The `alpha` value for trace (same across all 3 plots). Default is 1.
    IgnoreLabel : Boolean, optional
        Determines whether to omit the label of a particular trace from the 
        legend. Default is False.
    mlabel : string or None
        Stands for 'manual label'; can be used to create a label for a given 
        trace that is independent of the simulation or time parameters.
    color : string representing a color
        If not None, this setting will override the color set by the color cycle.
    markers : Boolean
        If True, the trace is plotted as markers instead of a continuous line.
    
    Returns
    -------
    A figure containing a set of plots represented on a linear, log-linear,
    and log-log scale.
    """
    
    if sys_obj is None:
        system_params = {}
    else:
        system_params = sys_obj.params().copy()
    
    if t_dict is None:
        t_dict = linear()
    else: 
        t_dict = t_dict.copy()
        
    time_params = {}
    for key, val in t_dict.items():
        if key == 'reprate' or key == 'N' or key == 'subsample' or key=='numcycles':
            time_params[key] = val
    
    varname = ivkey

    # preparing strings that will appear on figure legend and text
    if ivtype == 'tm':
        var = time_params.pop(ivkey)

        if varname == 'reprate':
            varname = 'reprate (MHz)'
            var_label = '%i' % var/1e6

        if varname == 'N' and var >= 1000:
            var_label = '%0.1fK' % (var/1000)
        
        
        else:
            var_label = '%i' % var

        rest_time_params = time_params
        rest_system_params = system_params

    elif ivtype == 'eq':
        var = system_params.pop(ivkey)
        var_label = '%0.2e' % var

        rest_time_params = time_params
        rest_system_params = system_params

    elif ivtype == 'other':
        var = other_params[varname]
        
        if varname == 'power':
            varname = 'Power (μW)'
            var_label = '%0.3f' % (var / 1e-6)
        
        else:
            var_label = '{}'.format(var)

        rest_time_params = time_params
        rest_system_params = system_params
    
    elif ivtype == 'none':
        varname = None
        var_label = None
        
        rest_time_params = time_params
        rest_system_params = system_params
    
    else:
        print("Invalid parameter type; must be 'tm', 'eq', 'other', or 'none' string.")
        return 0
    
    #--- Specifying unit
    unit_name = unit
    unit_value = units[unit]
    
    if annotate:
        txt = []
        for key, val in rest_time_params.items():
            if key == 'reprate':
                val /= 1e6
                p = 'reprate (MHz) = {0}'.format(int(val))

            elif key == 'N':
                val /= 1000
                val = int(val)
                p = ' {0} = {1}K'.format(key, val)
                
            else:
                p = '{0} = {1}'.format(key, val)
            
            txt.append(p)
        
        nL = True  # new line that separates time_params from system_params display
        
        for key, val in rest_system_params.items():
            p = '{0} = {1:0.2e}'.format(key, val)
            if nL:
                txt.append('\n'+p)
                nL = False
            else:
                txt.append(p)
        
        #print(txt)
    
    # Figure parameters
    fw = 13; fh = 5
    
    if fig is None:
        fresh_fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(fw, fh))

    else:
        fresh_fig = fig
        (ax1, ax2, ax3) = fig.axes
    
    if offset is None:
        offset = [0, 0 , 0 ]

    for i, ax in enumerate([ax1, ax2, ax3]):
        
        if ResetColorCyc:
            ax.set_prop_cycle(None)
            
        if mlabel is None:
            label = var_label
            leg_title = varname
        
        else:
            label = mlabel
            leg_title = None
            
        if IgnoreLabel:
            label = None
        
        xoff = x + offset[i]
        xoff /= unit_value

        if not norm:
            if markers:
                ax.plot(xoff, y, 'o', ms=4, 
                        mfc=ax.get_facecolor(), 
                        color=color, label=label, 
                        alpha=opaq)
            else:
                ax.plot(xoff, y, color=color, label=label, alpha=opaq,
                        linewidth=linewidth)

        else:
            if markers:
                ax.plot(xoff, kin_kit.normalized(y), 'o', ms=4, 
                        mfc=ax.get_facecolor(),
                        color=color, label=label, alpha=opaq)
            else:
                ax.plot(xoff, kin_kit.normalized(y), 
                    color=color, label=label, alpha=opaq, linewidth=linewidth)

        if fig is None:
            ax.set_xlabel('Time (%s)'%unit_name)
            if norm:
                ax.set_ylabel('Rel. Intensity')
            else:
                ax.set_ylabel('Counts')
                
        leg = ax1.legend(title=leg_title, fontsize=8, loc=1)
        plt.setp(leg.get_title(), fontsize=10)
    
    if xmin is None:
        xmin = 0.1
    if xmax is None:
        xmax = max(x)/unit_value
    
    ax1.set_xlim(0, xmax)
    ax2.set_xlim(0, xmax) 
    ax3.set_xlim(xmin, xmax)
    
    if norm:
        if ymax is None:
            ymax = 1.20
        if ymin is None:
            ymin = 0.001
            
        ax1.set_ylim(-0.05, 1.05)
        ax2.set_ylim(ymin, ymax) #0.001
        ax3.set_ylim(ymin, ymax)
    
    else:
        # set default ymax, ymin according to data limits
        if ymax is None: 
            # check if previously plotted trace has a higher maximum; 
            # if not, replace the current ymax with max(y)
            ymax = max([max(y), find_max_value(ax2)[1]])
        if ymin is None:                
            # analogous check for ymin
            ymin = min([min(y), find_min_value(ax2)[1]])
            
        ax1.set_ylim(-0.05*(ymax-ymin), 1.05*ymax)
        ax2.set_ylim(0.1*ymin, 2*ymax)
        ax3.set_ylim(0.1*ymin, 2*ymax)
    
    ax2.set_yscale('log')
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    
    if fig is None:
        fresh_fig.tight_layout()
    
    if annotate:
        fresh_fig.text(0, -0.2, ' , '.join(txt), transform = ax1.transAxes,
                 fontsize=8, color='darkgrey', fontstyle='italic')
        fresh_fig.subplots_adjust(bottom=0.202)
    
    return fresh_fig
        
def show_species(time, species_sets, system, norm, save=False, filename=None, destfolder=None):
    """
    Application of the plot3scales.plot function targeted at plotting species
    populations
    
    """
    species_names = system.populations
    alphas = np.linspace(0.3, 1, num=len(species_sets))
    if len(alphas)==1:
        alphas[0] = 1
    fig=None
    for k, species_set in enumerate(species_sets):
        for i, species in enumerate(species_set):
            species = kin_kit.align_by_max(species, time, species_sets[0][0], avgnum=5, value=0.5*units['ns'])
            fig=plot(time, species, sys_obj=None, t_dict=None,
                     annotate=False, fig=fig, linewidth=2, 
                     ResetColorCyc=i==0, opaq = alphas[k], norm=norm,
                     mlabel=species_names[i])
    plt.tight_layout()
    plt.show()
    
    if save:
        kin_kit.save_species(filename, destfolder, norm=norm)
        
def find_min_value(axis, which='y'):
    """
    Finds and returns the minimum y value among all plotted traces in an axis
    object.
    
    
    Parameters:
    -----------
    axis : axis object
        The axis (plot) of which we are trying to obtain the minimum value
    which : string, 'x' or 'y', optional
        Whether to look for the minimum x or y value. Default is y
    
    Returns:
    --------
    datapoint : tuple
        The x and y value of the found datapoint. Warning: if the minimum y 
        value is sought but there are multiple occurences of the same minimum 
        value, this function will return an x,y tuple where the x value 
        corresponds to the FIRST occurence of the minimum y value.
    
    """
    global_min = np.inf
    global_idx = None
    global_line = None
    
    lines = axis.lines
    
    if len(lines) == 0:
        return (np.inf, np.inf)
    
    for line in lines:
        # Step 1: find the local minimum of each line
        if which == 'y':
            local_min = line.get_ydata().min()
            local_idx = np.argmin(line.get_ydata())
        else:
            local_min = line.get_xdata().min()
            local_idx = np.argmin(line.get_xdata())
        # Step 2: check it against the current global minimum
        if local_min < global_min:
            global_min = local_min
            global_idx = local_idx
            global_line = line
        # Step 3: repeat for all lines
    
    # Step 4: Generate tuple
    if which == 'y':
        x_val = global_line.get_xdata()[global_idx]
        y_val = global_min
    else:
        x_val = global_min
        y_val = global_line.get_ydata()[global_idx]
    return x_val, y_val

def find_max_value(axis, which='y'):
    """
    Exact analogue of ``find_min_value()`` function, but for obtaining the 
    global maximum of the plots in a given axis object.
    
    """
    global_max = -np.inf
    global_idx = None
    global_line = None
    
    lines = axis.lines
    
    if len(lines) == 0:
        return (-np.inf, -np.inf)
    
    for line in lines:
        # Step 1: find the local minimum of each line
        if which == 'y':
            local_max = line.get_ydata().max()
            local_idx = np.argmax(line.get_ydata())
        else:
            local_max = line.get_xdata().max()
            local_idx = np.argmax(line.get_xdata())
        # Step 2: check it against the current global minimum
        if local_max > global_max:
            global_max = local_max
            global_idx = local_idx
            global_line = line
        # Step 3: repeat for all lines
    
    # Step 4: Generate tuple
    if which == 'y':
        x_val = global_line.get_xdata()[global_idx]
        y_val = global_max
    else:
        x_val = global_max
        y_val = global_line.get_ydata()[global_idx]
    return x_val, y_val
    
