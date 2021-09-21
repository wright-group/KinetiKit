KinetiKit
---------

INSTALLATION INSTRUCTIONS:

The latest releases of this package can be found on PyPI: https://pypi.org/project/KinetiKit/. 
An Anaconda Python installation is recommended: https://www.anaconda.com/distribution/
Open the Anaconda Command Prompt and run ``pip install KinetiKit``. To update to the latest version if you have already installed the package, run ``pip install KinetiKit -upgrade`` instead.

After installing the package, the simplest way to get started with the program is to come back to this Github page and download some of the files in the ``examples`` folder. By saving one of those files as a .py or .ipynb file in a directory that contains your experimental data, you can define file paths in the scripts relative to the path of the scripts themselves.

LIST OF EXAMPLE FILES

* calling_a_custom_model
* custom_models
* FitRates_Biexp
* FitRates_Hetero
* FitRates_Mono
* FitRates_Multipower
* FitRates_shorterTime
* imageShow
* plotTRPL3scales
* pulseViewer
* simplePlotPL
* simplePlotTRPL
* VisualizeSteadyState

In addition, a few iPython notebook files are available for interactive fitting:

* InteractivePlotting_Hetero
* InteractivePlotting_Mono
* InteractivePlotting_Multiexp

QUICK NOTES FOR NON-PYTHONIC FOLKS:

* Note on File paths: If ``dir_path = os.path.dirname(os.path.realpath(__file__))``, then all file paths will be read relative to where this file is located. Alternatively, you can set a directory path as it "normally" appears, e.g. ``dir_path= r'C:\Users\Me\Documents\MyData'``. In general, you can combine paths together (e.g. a path, a subfolder (or multiple subfolders) and a file name by using the command ``os.path.join(path1, path2, path3...)``.
* Note on Units: Many units and a few constants have been pre-defined in this package, inside ``KinetiKit.units``. On the provided example files, only the units that are used are imported. If you frequently use a unit that is not already imported, you can add it to the imported list at the top of the file (e.g. ``from units import ps, fs, MHz``).
* Note on Setting and Fitting parameters: As long as ``bounds`` have been set for a parameter, its value given in ``initparams`` does NOT affect the Differential Evolution search. Parameters for which ``bounds`` are not defined will be fixed at the value given in ``initparams``.
