#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 18:18:11 2023

@author: joe
"""

import sys
import pyscf
import multiprocessing
from functools import partial
import os
from pyscf.geomopt.geometric_solver import optimize
import time
from pyscf.hessian import thermo
import numpy as np


'''
Tests if a molecule has a singlet or a triplet ground state. This is done by 
calculating the energy of both singlet and triplet ground states and then seeing 
which is lower. 
Returns 1 if Singlet is lower, 3 if Triplet is lower, and 0 if there's an error'
'''
def check_multiplicity(path, name, basis_set, functional, verbose=False, initialchg=0):
    try:
        #Build singlet molecule
        mol_s = pyscf.M(
            atom = path+name,
            basis = basis_set,
            charge = initialchg,
            ecp = "cc-pvtz-pp", #Only used for heavy atoms
            spin = 0 #Same number of alpha and beta electrons
            )
        #Calculate singlet energy
        sing = pyscf.dft.UKS(mol_s)
        sing.xc = functional
        energy_s = sing.kernel()
    except Exception as e:
        print('Singlet calculation failed for ' + name + ':')
        print(e)
        return 0
    try:
        #Build triplet molecule
        mol_t = pyscf.M(
            atom = path+name,
            basis = basis_set,
            charge = initialchg,
            ecp = "cc-pvtz-pp", #Only used for heavy atoms
            spin = 2 #Two more alpha electrons than beta electrons
            )
        #Calculate triplet energy
        trip = pyscf.dft.UKS(mol_t)
        trip.xc = functional
        energy_t = trip.kernel()
    except Exception as e:
        print('Triplet calculation failed for ' + name + ':')
        print(e)
        return 0
    if energy_s < energy_t:
        return 1
    else:
        return 3

def make_mf(mol, functional):
    if functional == 'hf':
        if mol.spin == 0:
            mf=pyscf.scf.RHF(mol)
        else:
            mf=pyscf.scf.UHF(mol)
    else:
        if mol.spin == 0:
            mf=pyscf.dft.RKS(mol)
        else:
            mf=pyscf.dft.UKS(mol)
    return mf


'''
Returns True if the associated molecule has imaginary frequencies and False
otherwise
'''
def has_im_freqs(mol, functional):
    mf = make_mf(mol, functional)
    energy_1 = mf.kernel()
    hess = mf.Hessian().kernel()
    freqs = thermo.harmonic_analysis(mf.mol, hess, imaginary_freq=True)
    for (freq, eigenvector) in zip(freqs['freq_wavenumber'], freqs['norm_mode']):
        freqconj = np.conjugate(freq)
        if freq != freqconj:
            eigenvector = eigenvector*0.02
            for k in range(0, len(mol.atom)):
                for i in range(0, 3):
                    mol.atom[k][1][i] += eigenvector[k][i]*0.02
            #Recalculate energy
            mf2 = make_mf(mol, functional)
            energy_2 = mf2.kernel()
            if energy_2 < energy_1:
                print(freq)
                return True
    return False



'''
Performs local minimization over atomic coordinates for some xyz file, then 
makes a new xyz file with the optimized atomic coordinates. Returns the optimized coordinates.
bfcs is (basis set, functional, charge, spin). 
NOTE: only works for neutral molecules
'''
def optimize_coords(path, name, optpath, optname_modification, bfcs):
    filename, extension = name.rsplit('.',1)
    optname = f"{filename}{optname_modification}.{extension}"
    #construct the molecule
    mol = pyscf.M(
        atom = path+name,
        basis = bfcs[0],
        charge = bfcs[2],
        ecp = "cc-pvtz-pp",
        spin = bfcs[3]
        )
    try:
        # If the spin is 0, use RKS to save time. Otherwise, use UKS.
        if bfcs[3] != 0:
            neu = pyscf.dft.UKS(mol)
        else:
            neu = pyscf.dft.RKS(mol)
        # Optimize the coordinates
        neu.xc = bfcs[1]
        optneu = optimize(neu, maxsteps=400)
        # Print the coordinates
        with open(optpath+optname, 'w+') as out:
            out.write(str(len(optneu._atom)) + '\n')
            out.write(optname + '\n')
            for a in optneu._atom:
                out.write(a[0] + ' ' + str(round(a[1][0]*pyscf.lib.param.BOHR,6)) + ' ' + str(round(a[1][1]*pyscf.lib.param.BOHR,6)) + ' ' + str(round(a[1][2]*pyscf.lib.param.BOHR,6)) + '\n')
        return (name, optneu._atom, True)
    except Exception as e:
        print('-------------')
        print(e)
        print('-------------')
        return (name, e, False)


'''
Calculated the adiabatic detachment energy from an optimized atomic structure.
bfcs is (basis set, functional, charge, spin). 
Returns (name, ADE (in eV), <multiplicity of final state OR 0 if calculation failed> )
'''
def calc_ade(path, name, bfcs):
    mol = pyscf.M(
        atom = path+name,
        basis = bfcs[0],
        charge = bfcs[2],
        ecp = "cc-pvtz-pp",
        spin = bfcs[3])
    try: #Calculate energy of original molecule/ion
        ori = pyscf.dft.UKS(mol)
        ori.xc = bfcs[1]
        energy_init = ori.kernel()
        if not ori.converged:
            print('Fixed point energy calculation not converged')
            return (name, 'no convergence', 0)
    except Exception as e:
        print('Fixed point energy failed for ' + name + ':')
        print(e)
        return (name, e, 0)
    try: #Remove an electron, and calculate the optimized energy
        upmultmol = pyscf.M( #the case where multiplicity increases by 1
            atom = path+name,
            basis = bfcs[0],
            charge = bfcs[2]+1,
            ecp = "cc-pvtz-pp",
            spin = bfcs[3]+1 #test the case where spin goes up by 1 first
            )
        upmult = pyscf.dft.UKS(upmultmol)
        upmult.xc = bfcs[1]
        optupmult = optimize(upmult, maxsteps=400)
        if has_im_freqs(optupmult, bfcs[1]):
            return (name + ' has imaginary frequencies (saddle point) at upmult', 0, 0)
        fixedup = pyscf.dft.UKS(optupmult)
        fixedup.xc = bfcs[1]
        energy_up = fixedup.kernel()
        if bfcs[3] > 0:
            downmultmol = pyscf.M( #the case where multiplicity decreases by 1
                atom = path+name,
                basis = bfcs[0],
                charge = bfcs[2]+1,
                ecp = "cc-pvtz-pp",
                spin = bfcs[3]-1 #test the case where spin goes down by 1
                )
            downmult = pyscf.dft.UKS(downmultmol)
            downmult.xc = bfcs[1]
            optdownmult = optimize(downmult, maxsteps=400)
            if has_im_freqs(optdownmult, bfcs[1]):
                return (name + ' has imaginary frequencies (saddle point) at downmult', 0, 0)
            fixeddown = pyscf.dft.UKS(optdownmult)
            fixeddown.xc = bfcs[1]
            energy_down = fixeddown.kernel()
            if not (fixeddown.converged or fixedup.converged):
                return (name + ' no final calculations converged.', 0, 0)
            elif not fixeddown.converged:
                return (name + ' no convergence in mult ' + str(bfcs[3]), 0, 0)
            elif not fixedup.converged:
                return (name + ' no convergence in mult ' + str(bfcs[3]+2), 0, 0)
            elif (energy_down - energy_init) < (energy_up - energy_init):
                return (name, (energy_down - energy_init)*27.211407, bfcs[3])
        if not fixedup.converged:
            return (name + ' no convergence in mult ' + str(bfcs[3]+2), 0, 0)
        else:
            return (name, (energy_up - energy_init)*27.211407, bfcs[3]+2)
    except Exception as e:
        print('-------------')
        print(e)
        print('-------------')
        return (name, e, 0)
      
    
'''
Just calculates the first vertical detachment energy. This is the difference in
energy between a structure and that same structure with one fewer electron. 
The spin of the final state can vary, so we just report back the lowest possible VDE.
'''
def calc_vde(path, name, bfcs):
    mol = pyscf.M(
        atom = path+name,
        basis = bfcs[0],
        charge = bfcs[2],
        ecp = "cc-pvtz-pp",
        spin = bfcs[3])
    try:
        ori = pyscf.dft.UKS(mol)
        ori.xc = bfcs[1]
        energy_init = ori.kernel()
        if not ori.converged:
            print('Fixed point energy calculation not converged')
            return (name, 'no convergence', 0)
        upmultmol = pyscf.M( #the case where multiplicity increases by 1
            atom = path+name,
            basis = bfcs[0],
            charge = bfcs[2]+1,
            ecp = "cc-pvtz-pp",
            spin = bfcs[3]+1 #test the case where spin goes up by 1 first
            )
        upmult = pyscf.dft.UKS(upmultmol)
        upmult.xc = bfcs[1]
        energy_up = upmult.kernel()
        if bfcs[3] > 0:
            downmultmol = pyscf.M( #the case where multiplicity decreases by 1
                atom = path+name,
                basis = bfcs[0],
                charge = bfcs[2]+1,
                ecp = "cc-pvtz-pp",
                spin = bfcs[3]-1 #test the case where spin goes down by 1
                )
            downmult = pyscf.dft.UKS(downmultmol)
            downmult.xc = bfcs[1]
            energy_down = downmult.kernel()
            if not (downmult.converged or upmult.converged):
                return (name + ' no final calculations converged.', 0, 0)
            elif not downmult.converged:
                return (name + ' no convergence in mult ' + str(bfcs[3]), 0, 0)
            elif not upmult.converged:
                return (name + ' no convergence in mult ' + str(bfcs[3]+2), 0, 0)
            elif (energy_down - energy_init) < (energy_up - energy_init):
                return (name, (energy_down - energy_init)*27.211407, bfcs[3])
        if not upmult.converged:
            return (name + ' no convergence in mult ' + str(bfcs[3]+2), 0, 0)
        else:
            return (name, (energy_up - energy_init)*27.211407, bfcs[3]+2)
    except Exception as e:
        print('-------------')
        print(e)
        print('-------------')
        return (name, e, 0)
        
'''
Calculates the electronic excitation energies using TDDFT. Returns the name of
the file, a list of electronic excitation energies, and then either 0 if the 
calculation did not converge or 1 if it did.
'''
def calc_excitations(path, name, bfcs, nstates):
    mol = pyscf.M(
        atom = path+name,
        basis = bfcs[0],
        charge = bfcs[2],
        ecp = "cc-pvtz-pp",
        spin = bfcs[3])
    try: #calcualte the first n excitations
        groundstate = pyscf.dft.UKS(mol)
        groundstate.xc = bfcs[1]
        groundstate.kernel()
        mytd = pyscf.tddft.TDDFT(groundstate)
        mytd.singlet = False
        mytd.nstate = nstates
        mytd.kernel()
        mytd.analyze()
        energies = mytd.e*27.211407
        if groundstate.converged:
            for k in range(nstates):
                if not mytd.converged[k]:
                    energies[k] = 0
            return (name, energies, 1)
        else:
            return (name, energies, 0)
    except Exception as e:
        print('-------------')
        print(e)
        print('-------------')
        return (name, e, 0)


def process_molecule(name: str, fn_name: str, path: str):
    return check_multiplicity(path, name, 'ccpvtz', fn_name)



if __name__ == "__main__":
    start_time = time.time()
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

