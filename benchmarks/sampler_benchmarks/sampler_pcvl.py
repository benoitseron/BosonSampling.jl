import quandelibc as qc
import perceval as pcvl
import numpy as np
import time
from scipy import diagonal, randn
from scipy.linalg import qr
import random

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

def Generating_Input(n ,m , modes = None):
    "This function randomly chooses an input with n photons in m modes."
    if modes == None :
        modes = sorted(random.sample(range(m),n))
    state = "|"
    for i in range(m):
        state = state + "0"*(1 - (i in modes)) +"1"*(i in modes)+ ","*(i < m-1)
    return pcvl.BasicState(state + ">")

time_cliffords = np.empty([n+1,1])
Sampling_Backend = pcvl.BackendFactory().get_backend("CliffordClifford2017")

for ind, reps in enumerate(nreps):
    matrices = []
    for i in range(reps):
        size = ind+1
        nth = 1
        matrices.append(pcvl.Matrix.random_unitary(size))

    input_state = Generating_Input(ind+1,ind+1)
    start_qc_cliffords = time.time()
    for matrix in matrices:
        Sampling_Backend(matrix).sample(input_state)
    end_qc_cliffords = time.time()

    time_cliffords[ind] = (end_qc_cliffords - start_qc_cliffords)/reps
    print(ind+1, time_cliffords[ind])

f = open("benchmarks/data/cliffords-pcvl.txt", "w")
for row in time_cliffords:
    np.savetxt(f, row)
f.close()
