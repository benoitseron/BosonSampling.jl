import thewalrus
from thewalrus import perm
import numpy as np
import time
from scipy import diagonal, randn
from scipy.linalg import qr

def haar_measure(n):
    '''A Random matrix distributed with Haar measure
    See https://arxiv.org/abs/math-ph/0609050
    How to generate random matrices from the classical compact groups
    by Francesco Mezzadri '''
    z = (np.random.randn(n,n) + 1j*np.random.randn(n,n))/np.sqrt(2.0)
    q,r = qr(z)
    d = np.diagonal(r)
    ph = d/np.abs(d)
    q = np.multiply(q,ph,q)
    return q

# benchmark inspired from https://the-walrus.readthedocs.io/en/latest/gallery/permanent_tutorial.html &
# https://github.com/Quandela/Perceval/blob/main/scripts/performance.py
a0 = 300.
anm1 = 2
# n = 29
n = 29
r = (anm1/a0)**(1./(n-1))
nreps = [(int)(a0*(r**((i)))) for i in range(n+1)]
times_walrus = np.empty([n+1,1])

for ind, reps in enumerate(nreps):
    matrices = []
    for i in range(reps):
        size = ind+1
        nth = 1
        matrices.append(haar_measure(size))
    start_walrus = time.time()
    for matrix in matrices:
        res = thewalrus.perm(matrix)
    end_walrus = time.time()

    times_walrus[ind] = (end_walrus - start_walrus)/reps
    print(ind+1, times_walrus[ind])

f = open("benchmarks/thewalrus_data.txt", "w")
for row in times_walrus:
    np.savetxt(f, row)
f.close()
