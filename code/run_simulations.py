import subprocess
import time
import sys

import numpy as np

from sympy.solvers import solve
from sympy import Symbol

from tqdm import tqdm


#NOTE: Probability and Costs parameters passed to simulation are interchanged, depending on specific model


def find_rhos(N_ICU, mu, rhos):
    x = Symbol('x')
    
    lrs = []
    for r in rhos:
        seqn = solve(x / (N_ICU * mu) - r, x)
        lrs.append((np.around(float(seqn[0]), decimals=3), np.around(r, decimals=1)))
        
    return lrs


def run_sim(ED, ICU, A, T, MU, icu_p, rr_p, beta, l, rho, counter, s_type='base'):
    if s_type == 'base':
        #./sim ED ICU A T lambda mu icup rrp rho counte
        subprocess.call('./base_sim ' + str(ED) + ' ' + str(ICU) + ' ' + str(A) + ' ' + str(T) + \
                        ' ' + str(l) + ' ' + str(MU) + ' ' + str(icu_p) + ' ' + str(rr_p) + ' ' + beta + ' ' + \
                        str(rho) + ' ' + str(counter), shell=True)
    elif s_type == 'dp':
        subprocess.call('./dp_sim ' + str(ED) + ' ' + str(ICU) + ' ' + str(A) + ' ' + str(T) + \
                        ' ' + str(l) + ' ' + str(MU) + ' ' + str(icu_p) + ' ' + str(rr_p) + ' ' + beta + ' ' + \
                        str(rho) + ' ' + str(counter), shell=True)
       

if __name__ == '__main__':
    ED = 15
    ICU = float(sys.argv[1])
    A = 3
    T = 8760
    MU = float(sys.argv[2])
    beta = sys.argv[3]
    phi = -3.0
    h = 3.0

    rhos = np.arange(0.1, 1.2, .1)
    lambda_rhos = find_rhos(ICU, MU, rhos)

    if sys.argv[4] == 'base':
        print(lambda_rhos)
        probs = [0.0, 0.25, 0.50, 0.75, 1.0]
        for tup in tqdm(lambda_rhos):
            for ip in probs:
                for rp in probs:
                    if ip == 1.0 and rp == 0.0:
                        counter = 0
                        while counter < 2:
                            #beta (discount factor) not needed in state independent model
                            run_sim(ED, ICU, A, T, MU, ip, rp, 0, tup[0], tup[1], counter, 'base')
                            counter += 1
    elif sys.argv[4] == 'dp':
        print(lambda_rhos)
        for tup in tqdm(lambda_rhos):
            counter = 0
            for counter in range(0, 10):
                run_sim(ED, ICU, A, T, MU, phi, h, beta, tup[0], tup[1], counter, 'dp')
                counter += 1
