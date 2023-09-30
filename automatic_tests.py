#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 27 14:57:27 2023

@author: joe
"""

import gases_functions
import os
import sys
import time

start_time = time.time()

# Define a few colors
reset = "\033[0m"
green = "\033[32m"
red = "\033[31m"

# Molecule files
O2_string = '''2
Oxygen
O -0.604112 0.0 0.0
O 0.604112 0.0 0.0
'''
H2O_string = '''3
Water
O 0.0 0.0 0.0
H 0.0 1.0 0.0
H 0.0 0.0 1.0
'''

# Define the directory for all molecules to be tested.
test_directory = 'test_molecules/'

# Define parameters for calculations (basis set, functional, initial charge)
basis_set = 'ccpvtz'
functional = 'b3lyp'
charge = 0

# Make the directory that holds test molecules
try:
    os.makedirs(test_directory)
except:
    print('Directory ' + test_directory + ' exists')
        

# Write the original files for the test molecules
with open(test_directory + 'O2.xyz', 'w+') as xyz:
    xyz.write(O2_string)

with open(test_directory + 'H2O.xyz', 'w+') as xyz:
    xyz.write(H2O_string)


# Optimize O2 and H2O's geometry
optname_modification = '_opt'
O2_opt = gases_functions.optimize_coords(test_directory, 'O2.xyz', test_directory, optname_modification, (basis_set, functional, charge, 2))
H2O_opt = gases_functions.optimize_coords(test_directory, 'H2O.xyz', test_directory, optname_modification, (basis_set, functional, charge, 0))

# Calculate VDEs
O2_vde = gases_functions.calc_vde(test_directory, 'O2'+optname_modification+'.xyz', (basis_set, functional, charge, 2))
H2O_vde = gases_functions.calc_vde(test_directory, 'H2O'+optname_modification+'.xyz', (basis_set, functional, charge, 0))

# Calculate ADEs
O2_ade = gases_functions.calc_ade(test_directory, 'O2'+optname_modification+'.xyz', (basis_set, functional, charge, 2))
H2O_ade = gases_functions.calc_ade(test_directory, 'H2O'+optname_modification+'.xyz', (basis_set, functional, charge, 0))

# Check whether O2 and H2O are calculated to have the correct multiplicities
O2_mult = gases_functions.check_multiplicity(test_directory, 'O2.xyz', basis_set, functional)

if O2_mult == 3:
    print(green + "Multiplicity of O2 matches!" + reset)
else:
    print(red + "O2 Multiplicity is incorrect!" + reset)
    
H2O_mult = gases_functions.check_multiplicity(test_directory, 'H2O.xyz', basis_set, functional)

if H2O_mult == 1:
    print(green + "Multiplicity of H2O matches!" + reset)
else:
    print(red + "H2O Multiplicity is incorrect!" + reset)

#Report values
print('Oxygen results:')
print('VDE: ')
print(O2_vde)
print('Experimental value is 12.3 eV')
print('Should get 12.86 eV for VDE')
print('b3lyp/ccpvtz benchmarked VDE is 12.96')
print('ADE')
print(O2_ade)
print('Should get 12.38 eV for ADE')
print('--------------------------------------------')
print('Water results:')
print('VDE: ')
print(H2O_vde)
print('Experimental value is 12.62 eV')
print('b3lyp/ccpvtz benchmarked VDE is 12.81')
print('ADE')
print(H2O_ade)

print("--- %s seconds ---" % (time.time() - start_time))











