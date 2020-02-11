KinetiKit

INSTALLATION INSRUCTIONS:
The latest releases of this package can be found on PyPI: https://pypi.org/project/KinetiKit/. 
An Anaconda Python installation is recommended: https://www.anaconda.com/distribution/
Open the Anaconda Command Prompt and run: 

pip install KinetiKit

To update to the latest version:

pip install KinetiKit -upgrade


LIST OF EXAMPLE FILES (INSTRUCTIONS BELOW)

* FitRates_Hetero
* FitRates_Mono
* plotTRPL3scales
* pulseViewer
* simplePlotPL
* simplePlotTRPL
* simRates_Hetero
* simRates_Mono
* VisualizeSteadyState

INSTRUCTIONS FOR FitRates_Hetero AND FitRates_Mono:
* Specify doFit; if False, the program will simply simulate the system with the current parameters.
	If True, the program will attempt a fit of simulation and data by varying parameters.
* Specify keepParams (usually False); if you set doFit = False and keepParams = True, the program
	will simulate the system with the parameters from the last time the program was run (e.g.
	from a previous fit), instead of with initparams.
* Define dir_path = r'C:\...' according to the directory of your data files
* Define subfolder as r'foldername', if it exists, otherwise set subfolder = ''
* Type the file name of your file(s) as a list of strings.
* If you don't want to follow the suggested file-name/key/title format, you can redefine the 
	keys list and the title string.
* In the line all_data.append(...) you may want to redefine t_zero, skip_h, and skip_f 
	according to the structure of your data file; t_zero is typically the time when counts 
	first get recorded, skip_h is the number of lines in the text file until that point, 
	and skip_f is the lines skipped from the bottom until the time point where a full 
	repetition cycle has occured (e.g. if your rep. rate is 80 MHz, then the final time
	point imported should be t_zero + 12.5 ns = 0.582259 + 12.5. 
* Define all initparams pertaining to your system of choice (e.g. for class Mono, the parameters
	are k_ann, k_dis, k_rec, and cs.
* Define the upper and lower bounds for only the parameters you want to be varied during the 
	simulation.
* Define the light and system instances; see sim.lib for default excitation parameters
	and sim.systems for options for systems. A "system" is essentially the model and its 
	corresponding rate equations that describe the behavior of charge carriers in the compound
	studied. A "light" object contains information about the power, repetition rate, wavelength,
	and pulse duration of the excitation light.
* In sim_comp_args, the relevant arguments you may want to change are:
	* irf_fwhm : the FWHM of your instrument response function
	* comparison : 'linear' or 'log' comparison between normalized data and simulation 
	* absolute : True or False; whether to consider the absolute difference between the data 
		and simulation arrays, or divide the difference by each data value (thus avoiding the
		over-weighing of small data points).
	* roll_value, roll_criterion : shifts ("rolls") both the data and simulation arrays so that
		their maxima (or steepest points, depending on roll_criterion) lie along roll_value in the
		time axis. It does not make much difference to roll to a value different than zero, unless
		you are using 'limits' to selectively fit part of the data. In which case, it makes sense to
		set roll_value to the same value as you set it later (aligned_data = TRPL_kit.align...)
		so that you can visualize where to set the limits.
	* limits : use it to selectively fit only a portion of the data. Careful; the values set as limits
		within sim_comp_args may not be the same as the values you see when the graph is plotted. To 
		avoid confusion, always set roll_value to the same value as set when aligning the data and sims
		later in the code (aligned_data = ... kin_kit.align....value=...*ns)
	* maxavgnum : if roll_criterion = 'max', maxavgnum determines how many highest data points are averaged
		to determine the "maximum value" of the arrays to be rolled.
* art.plot3scales.plot follows a default color scheme and labeling system, but you can choose several parameters
	to overwrite the defaults:
	* Plotting parameters: linewidth, markers, color
	* mlabel : by default, each trace will be labeled by its key as defined above, and either '_data' or '_sim'.
* kin_kit.save_fit(...) : this line is commented out unless you want to save the traces and parameters of a
	recent fit. You can copy and paste this line directly to the iPython console and run it after the fit is
	finished, or un-comment it so that the files are automatically saved once the fitting is complete. Note 
	that this does NOT save the produced figure.

INSTRUCTIONS FOR simplePlotPL AND simplePlotTRPL:
These are simple example files simply for plotting data. You simply need to define the file path, subfolder,
filenames, and identifiers (keys) for each file name, as well as a title for the output plot.

Coming soon:
INSTRUCTIONS FOR simRates_Hetero AND simRates_Mono
INSTRUCTIONS FOR VisualizeSteadyState
INSTRUCTIONS FOR pulseViewer
