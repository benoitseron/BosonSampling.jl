import thewalrus
from thewalrus import perm
import perceval as pcvl
import quandelibc as qc
import numpy as np
import time
from scipy import diagonal, randn
from scipy.linalg import qr

def haar_measure(n):
    '''A Random matrix distributed with Haar measure
    See https://arxiv.org/abs/math-ph/0609050
    How to generate random matrices from the classical compact groups
    by Francesco Mezzadri '''
    z = (randn(n,n) + 1j*randn(n,n))/np.sqrt(2.0)
    q,r = qr(z)
    d = diagonal(r)
    ph = d/np.abs(d)
    q = np.multiply(q,ph,q)
    return q

# benchmark inspired from https://the-walrus.readthedocs.io/en/latest/gallery/permanent_tutorial.html &
# https://github.com/Quandela/Perceval/blob/main/scripts/performance.py
a0 = 300.
anm1 = 2
n = 30
r = (anm1/a0)**(1./(n-1))
nreps = [(int)(a0*(r**((i)))) for i in range(n)]
times_walrus = np.empty(n)
times_qc_1 = np.empty(n)

for ind, reps in enumerate(nreps):
    #print(ind+1,reps)
    matrices = []
    for i in range(reps):
        size = ind+1
        nth = 1
        # matrices.append(pcvl.Matrix.random_unitary(size))
        matrices.append(haar_measure(size))
    start_walrus = time.time()
    for matrix in matrices:
        res = thewalrus.perm(matrix)
    end_walrus = time.time()
    start_qc_1 = time.time()
    for matrix in matrices:
        res = qc.permanent_cx(matrix, 2)
    end_qc_1 = time.time()

    times_walrus[ind] = (end_walrus - start_walrus)/reps
    times_qc_1[ind] = (end_qc_1 - start_qc_1)/reps

    print(ind+1, times_walrus[ind], times_qc_1[ind])
    # print(ind+1, times_walrus[ind])

res_thewalrus = list(times_walrus)
with open("data_thewalrus.txt", 'w') as f:
    for res in res_thewalrus:
        f.write("%s " %res)
res_pcvl = list(times_qc_1)
with open("data_pcvl.txt", 'w') as f:
    for res in res_pcvl:
        f.write("%s " %res)
