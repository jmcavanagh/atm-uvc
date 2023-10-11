
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  6 17:53:52 2023

@author: joe
"""

import gases_functions
import multiprocessing
import sys
import os
from functools import partial
import time
import pyscf

def check_multiplicity_anion(path, name, basis_set, functional, verbose=False, initialchg=1):
    try:
        #Build doublet molecule
        mol_d = pyscf.M(
            atom = path+name,
            basis = basis_set,
            charge = initialchg,
            ecp = "cc-pvtz-pp", #Only used for heavy atoms
            spin = 1 #One more alpha than beta electron
            )
        #Calculate doublet energy
        doub = pyscf.dft.UKS(mol_d)
        doub.xc = functional
        energy_d = doub.kernel()
    except Exception as e:
        print('Doublet calculation failed for ' + name + ':')
        print(e)
        return 0
    try:
        #Build quartet molecule
        mol_q = pyscf.M(
            atom = path+name,
            basis = basis_set,
            charge = initialchg,
            ecp = "cc-pvtz-pp", #Only used for heavy atoms
            spin = 3 #Three more alpha electrons than beta electrons
            )
        #Calculate triplet energy
        quar = pyscf.dft.UKS(mol_q)
        quar.xc = functional
        energy_q = quar.kernel()
    except Exception as e:
        print('Quartet calculation failed for ' + name + ':')
        print(e)
        return 0
    if energy_d < energy_q:
        return 2
    else:
        return 4



def process_molecule(name, path, bfcs):
    multiplicity = check_multiplicity_anion(path, name, bfcs[0], bfcs[1])
    print(name + ':')
    print(multiplicity)
    if multiplicity == 2:
        bfcs[3] = 1
    elif multiplicity == 4:
        bfcs[3] = 3
    else:
        return 'something went wrong. Multiplicity = ' + str(multiplicity)
    mol = pyscf.M(
        atom = path+name,
        basis = bfcs[0],
        charge = bfcs[2],
        ecp = "cc-pvtz-pp",
        spin = bfcs[3]
        )
    imfreqs = gases_functions.has_im_freqs(mol, bfcs[1])
    if imfreqs == True:
        print(path+name + ' has imaginary frequencies!')
        return path+name
    else:
        print(path+name + ' has real frequencies')
        return 'real'


if __name__ == "__main__":
    folder = sys.argv[1]

    files = os.listdir(folder)


    start_time = time.time()
    pool = multiprocessing.Pool(processes=len(files))
    plus_args = partial(process_molecule, path=folder, bfcs = ['sto3g', 'cam-b3lyp', 1, 1])

    results = pool.map(plus_args, files)

    for r in results:
        if r != 'real':
            print(r)
    print('------------')
    print("--- %s seconds ---" % (time.time() - start_time))
