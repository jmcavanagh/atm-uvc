import gases_functions
import pyscf
import multiprocessing
import sys
import os
from functools import partial
import time
import numpy as np


# Returns the largest force felt by a single atom in units of Hartree/Bohr
def get_largest_force(name, path, bfcs):
    multiplicity = gases_functions.check_multiplicity(path, name, bfcs[0], bfcs[1])
    bfcs[3] = multiplicity - 1
    mol = pyscf.M(
        atom = path+name,
        basis = bfcs[0],
        symmetry = False,
        charge=bfcs[2],
        spin=bfcs[3])
    g = pyscf.dft.UKS(mol).run(xc=bfcs[1]).nuc_grad_method()
    grads = g.kernel()
    return np.max(np.sum(grads**2,1)**0.5)


if __name__ == "__main__":
    folder = sys.argv[1]

    files = os.listdir(folder)

    fn_names = ['b3lyp']

    for fn_name in fn_names:
        start_time = time.time()
        pool = multiprocessing.Pool(processes=len(files))
        plus_args = partial(get_largest_force, path=folder, bfcs = ['ccpvtz', fn_name, 0, 0])

        results = pool.map(plus_args, files)

        print(fn_name)
        for r in results:
            print(r)
        print('------------')
        print("--- %s seconds ---" % (time.time() - start_time))
