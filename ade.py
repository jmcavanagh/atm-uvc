#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 18:18:11 2023

@author: joe
"""

import sys
import pyscf
from pyscf import gto, scf, dft, tddft, cc
from multiprocessing import Process, Queue
import multiprocessing
from functools import partial
import os
from pyscf.geomopt.geometric_solver import optimize
import time
start_time = time.time()

def get_file_coords(path, name):
    with open(path+name) as fil:
        molec_structure = ''
        for line in fil:
            if not line[0].isdigit():
                molec_structure += line
    return molec_structure


def calculate_ie_dft(path, name, basis_set, functional, verbose=False, initialchg=0):
    molec_structure = get_file_coords(path, name)
    #construct the molecule
    mol = pyscf.M(
        atom = molec_structure,
        basis = basis_set,
        charge = initialchg
        )
    try:
        # Optimize the neutral structure
        neu = pyscf.dft.RKS(mol)
        neu.xc = functional
        optneu = optimize(neu, maxsteps=400)
        # Now,set up a fine calculation of energy
        nf = pyscf.dft.UKS(optneu)
        energy_i = nf.kernel()
        # Now, remove an electron.
        optneu.charge += 1
        optneu.spin = 1
        cat = pyscf.dft.UKS(optneu)
        # Now, optimize the cation structure
        optcat = optimize(cat, maxsteps=400)
        # Do another fine calculation of energy
        cf = pyscf.dft.UKS(optcat)
        energy_f = cf.kernel()
        # Return difference between optimized molecule and optimized cation, in eV
        E = (energy_f - energy_i)*27.211407
        # Check that the calculations converged and return
        if cf.converged and nf.converged:
            return (name, E, True)
        return (name, E, False)
    except Exception as e:
        print('-------------')
        print(e)
        print('-------------')
        return (name, 0, False)


def process_molecule(name: str, fn_name: str, path: str):
    return calculate_ie_dft(path, name, 'ccpvtz', fn_name)


folder = sys.argv[1]

files = os.listdir(folder)

fn_names = ['b3lyp', 'pbe0', 'cam-b3lyp']


for fn_name in fn_names:
    pool = multiprocessing.Pool(processes=len(files))

    plus_args = partial(process_molecule, fn_name = fn_name, path=folder)

    results = pool.map(plus_args, files)

    print(fn_name)
    for r in results:
        print(r)
    print('------------')
    print("--- %s seconds ---" % (time.time() - start_time))
