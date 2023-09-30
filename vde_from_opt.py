#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 00:44:02 2023

@author: joe
"""


import gases_functions
import multiprocessing
import sys
import os
from functools import partial
import time


def process_molecule(name, path, bfcs, logpath):
    logfile = name + '_' + bfcs[1] + '_vde.log'
    sys.stdout = open(logpath+logfile, 'w+')
    multiplicity = gases_functions.check_multiplicity(path, name, bfcs[0], bfcs[1])
    print(name + ':')
    print(multiplicity)
    print('------')
    if multiplicity == 1:
        bfcs[3] = 0
    elif multiplicity == 3:
        bfcs[3] = 2
    else:
        sys.stdout.close()
        return 'something went wrong. Multiplicity = ' + multiplicity
    vde_tuple = gases_functions.calc_vde(path, name, bfcs)
    sys.stdout.close()
    return vde_tuple

if __name__ == "__main__":
    folder = sys.argv[1]

    files = os.listdir(folder)

    #fn_names = ['b3lyp', 'pbe0', 'cam-b3lyp']
    fn_names = ['b3lyp']

    for fn_name in fn_names:
        start_time = time.time()
        pool = multiprocessing.Pool(processes=len(files))
        plus_args = partial(process_molecule, path=folder, bfcs = ['ccpvtz', fn_name, 0, 0], logpath = './logs/')

        results = pool.map(plus_args, files)

        print(fn_name)
        for r in results:
            print(r)
        print('------------')
        print("--- %s seconds ---" % (time.time() - start_time))
