import strawberryfields as sf
from strawberryfields.ops import *
import numpy as np
from numpy import pi, sqrt
import scipy as sp
from scipy.stats import unitary_group
import scipy.linalg as la
import matplotlib.pyplot as plt
import timeit
from tqdm import tqdm

def QFT(n):
    F = np.zeros((n,n), dtype=complex)
    for i in range(n):
        for j in range(n):
            F[i,j] = np.exp(-2j*np.pi/n)**(i*j)
    return 1/(np.sqrt(n)) * F

def random_unitary(n):
    U = unitary_group.rvs(n)
    return np.dot(U, U.conj().T)

def sampler(n):
    boson_sampling = sf.Program(n)
    U = QFT(n)

    with boson_sampling.context as q:
        for i in range(n):
            Fock(1) | q[i]
        Interferometer(U) | q
        # MeasureFock() | q

    eng = sf.Engine(backend="fock", backend_options={"cutoff_dim": n})
    eng.run(boson_sampling)

res_python = []
for n in tqdm(range(2,10+1), ncols=70):
    res = timeit.timeit('sampler(n)', globals=globals(), number=1)
    res_python.append(res)

with open("data_python.txt", 'w') as f:
    for res in res_python:
        f.write("%s " %res)
