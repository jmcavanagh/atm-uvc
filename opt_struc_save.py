
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




def process_molecule(name, path, optname_modification, optpath, bfcs, logpath):
    logfile = name + '_' + bfcs[1] + '_opt.log'
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
        return 'something went wrong. Multiplicity = ' + multiplicity
    opt_tuple = gases_functions.optimize_coords(path, name, optpath, optname_modification, bfcs)
    sys.stdout.close()
    return opt_tuple


if __name__ == "__main__":
    folder = sys.argv[1]

    files = os.listdir(folder)

    #fn_names = ['b3lyp', 'pbe0', 'cam-b3lyp']
    #fn_names = ['b3lyp']
    fn_names = ['pbe0', 'cam-b3lyp']
    
    for fn_name in fn_names:
        start_time = time.time()
        pool = multiprocessing.Pool(processes=len(files))
        plus_args = partial(process_molecule, path=folder, optname_modification='_opt'+fn_name, optpath='./optstrucs/', bfcs = ['ccpvtz', fn_name, 0, 0], logpath='./logs/')

        results = pool.map(plus_args, files)

        print(fn_name)
        for r in results:
            print(r)
        print('------------')
        print("--- %s seconds ---" % (time.time() - start_time))

