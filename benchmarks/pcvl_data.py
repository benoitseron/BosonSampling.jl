import quandelibc as qc
import perceval as pcvl
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
n = 29
r = (anm1/a0)**(1./(n-1))
nreps = [(int)(a0*(r**((i)))) for i in range(n+1)]
times_qc_glynn = np.empty([n+1,1])
times_qc_ryser = np.empty([n+1,1])

for ind, reps in enumerate(nreps):
    matrices = []
    for i in range(reps):
        size = ind+1
        nth = 1
        matrices.append(haar_measure(size))

    start_qc_glynn = time.time()
    for matrix in matrices:
        res = qc.permanent_cx(matrix, 2)
    end_qc_glynn = time.time()

    start_qc_ryser_4 = time.time()
    for matrix in matrices:
        res = qc.permanent_cx(matrix, 4)
    end_qc_ryser_4 = time.time()

    times_qc_glynn[ind] = (end_qc_glynn - start_qc_glynn)/reps
    times_qc_ryser[ind] = (end_qc_ryser_4 - start_qc_ryser_4)/reps
    print(ind+1, times_qc_glynn[ind], times_qc_ryser)

f = open("benchmarks/glynn-pcvl.txt", "w")
for row in times_qc_glynn:
    np.savetxt(f, row)
f.close()

f = open("benchmarks/ryser4-pcvl.txt", "w")
for row in times_qc_ryser:
    np.savetxt(f, row)
f.close()
