# atm-uvc
Studying the effects of Far-UVC radiation on atmospheric components

This gives a toolbox to calculate the adiabatic detachment energies, vertical detachment energies, and electronic excitation energies of a wide range of molecules. Primarily designed for gases, but can work for a wide range of other molecules.

The primary file that contains functions used for all other files is `gases_functions.py`.

To use this, start with a directory (or several) of .xyz files for all molecules you're interested in.

To optimize structures in a folder:

Open `opt_struc_save.py` and edit the variable named `fn_names` to include all exchange-correlation functionals you want to use.
Make a directory named "logs":
`mkdir logs`
Make a directory named "optstrucs" (or edit the variable optpath in the python file to another directory to store optimized structures).
`mkdir optstrucs`
Then, run `python3 opt_struc_save.py <folder>` to optimize all of the structures in `<folder>` and save the optimized structures to the folder `optstrucs`.

To calculate adiabatic detachment energy:

Make a directory named "logs":
`mkdir logs`



