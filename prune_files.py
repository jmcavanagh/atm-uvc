#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 01:28:20 2023

@author: joe
"""

import multiprocessing
import sys
import os
from functools import partial
import time
import pyscf
import numpy as np

def atoms_too_close(testmol, name):
    try:
        testmol.build()
        coords = testmol.atom_coords(unit='Bohr')
        for k in range(0, len(coords)):
            for l in range(k+1, len(coords)):
                if np.sum((coords[k] - coords[l])**2) < 0.99:
                    print('---')
                    print(coords[k])
                    print(coords[l])
                    print(np.sum((coords[k] - coords[l])**2))
                    print('---')
                    return True
    except Exception as e:
        print(name)
        print(e)
        return True
    return False

def process_molecule(name, path):
    try:
        testmol = pyscf.gto.Mole(atom = path+name)
        if testmol.nelectron % 2 != 0:
            print(name + " has odd number of electrons")
            return name
        elif atoms_too_close(testmol, name):
            print(name + " has atoms too close together")
            return name
        else:
            return 0
    except Exception as e:
        print(name)
        print(e)
        return name
    
    
if __name__ == "__main__":
    num_bad_files = 0
    num_files = 0
    rmcommand = 'rm '
    for n in range(1, 31):
        folder = '../volatile_hs_pruned_splitup/folder_' + str(n) + '/'
        
        files = os.listdir(folder)
        
        pool = multiprocessing.Pool(processes=len(files))
        plus_args = partial(process_molecule, path=folder)
    
        results = pool.map(plus_args, files)
        for r in results:
            if r != 0:
                rmcommand += r + ' '
                num_bad_files += 1
            num_files += 1
        
    print(rmcommand)
    print(num_bad_files)
    print(num_files)
    print(num_bad_files / num_files)










