
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





def process_molecule(name, path, bfcs):
    multiplicity = gases_functions.check_multiplicity(path, name, bfcs[0], bfcs[1])
    print(name + ':')
    print(multiplicity)
    if multiplicity == 1:
        bfcs[3] = 0
    elif multiplicity == 3:
        bfcs[3] = 2
    else:
        return 'something went wrong. Multiplicity = ' + multiplicity
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

    #fn_names = ['b3lyp', 'pbe0', 'cam-b3lyp']
    fn_names = ['b3lyp']
    #fn_names = ['pbe0', 'cam-b3lyp']
    
    for fn_name in fn_names:
        start_time = time.time()
        pool = multiprocessing.Pool(processes=len(files))
        plus_args = partial(process_molecule, path=folder, bfcs = ['ccpvtz', fn_name, 0, 0])

        results = pool.map(plus_args, files)

        print(fn_name)
        for r in results:
            if r != 'real':
                print(r)
        print('------------')
        print("--- %s seconds ---" % (time.time() - start_time))

