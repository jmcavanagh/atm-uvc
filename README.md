# atm-uvc
Studying the effects of Far-UVC radiation on atmospheric components

This gives a toolbox to calculate the adiabatic detachment energies, vertical detachment energies, and electronic excitation energies of a wide range of molecules. Primarily designed for gases, but can work for a wide range of other molecules.

#Theory
##On ionization
Removing an electron from a molecule requires some amount of energy, which is called the ionization energy. During the process of removing an electron, the structure of a molecule can change. After the electron is removed, the final structure of the molecule can very, and so the energy cost of the electron removal can also change. The most likely ionization energy, called the 'vertical detachment energy' or 'VDE' involves no change between the structures. The VDE of a molecule is simply calculated as the difference in energy between a molecule with some structure and an ion with the same structure as the molecule, but missing an electron. The lowest possible ionization energy is called "adiabatic detachment energy" or "ADE" and is the difference between the molecule's energy and the lowest possible ion's energy.

All ADE's and VDE's calculated by this program use density functional theory.

##On excitation
Absorption of a photon can change the energy of a molecule by exciting an electron to a different orbital. The energies possible for this can be relatively well-approximated using time-dependent density functional theory (TDDFT)

##On imaginary frequencies and saddle points
The energy of some molecule depends on the positions of all of its atomic nuclei. Therefore, the energy can be thought of as a function over all of the positions of all of the nuclei. Local minima of energy represent stable molecular structures, so we can often generate a more accurate structure from some template using local minimization methods. However, some molecular structures can be at saddle points, and these methods will fail. We can check if a molecule is at a saddle point by seeing if its hessian has negative eigenvalues, which correspond to imaginary vibrational frequencies. The files `check_for_im_freqs.py` and `check_for_im_freqs_anions.py` serve the purpose of searching for saddle points.


#How to use

The primary file that contains functions used for all other files is `gases_functions.py`.

To use this, start with a directory (or several) of .xyz files for all molecules you're interested in.

##To optimize structures in a folder:

Open `opt_struc_save.py` and edit the variable named `fn_names` to include all exchange-correlation functionals you want to use.

Make a directory named "logs":
`mkdir logs`

Make a directory named "optstrucs" (or edit the variable optpath in the python file to another directory to store optimized structures).
`mkdir optstrucs`

Then, run `python3 opt_struc_save.py <folder>` to optimize all of the structures in `<folder>` and save the optimized structures to the folder `optstrucs`.

##To calculate adiabatic detachment energy:

Open `ade_from_opt.py` and edit the variable named `fn_names` to include all exchange-correlation functionals you want to use.

Make a directory named "logs":
`mkdir logs`

Then, run `python3 ade_from_opt.py <folder>` to optimize all of the structures in `<folder>` and save the optimized structures to the folder `optstrucs`.

Generally, you should optimize both the neutral and cation structures with the same exchange-correlation functional. What this means is that you should only calculate the ADE using "b3lyp" if the starting structure was optimized with "b3lyp".

## To check

