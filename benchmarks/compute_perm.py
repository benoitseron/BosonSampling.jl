import csv
from scipy.stats import unitary_group
import numpy as np
import timeit
from tqdm import tqdm

def ryser(M):

    n = len(M)
    row_comb = np.zeros((n), dtype=M.dtype)

    total = 0
    old_grey = 0
    sign = +1
    binary_power_dict = [2**i for i in range(n)]
    num_loops = 2**n

    for k in range(0, num_loops):
        bin_index = (k + 1) % num_loops
        reduced = np.prod(row_comb)
        total += sign * reduced
        new_grey = bin_index ^ (bin_index // 2)
        grey_diff = old_grey ^ new_grey
        grey_diff_index = binary_power_dict.index(grey_diff)
        new_vector = M[grey_diff_index]
        direction = (old_grey > new_grey) - (old_grey < new_grey)
        for i in range(n):
            row_comb[i] += new_vector[i] * direction
        sign = -sign
        old_grey = new_grey

    return total

res_python = []
for n in tqdm(range(15,25+1), ncols=70):
    U = unitary_group.rvs(n)
    res = timeit.timeit('ryser(U)', globals=globals(), number=1)
    res_python.append(res)

with open("data_python.txt", 'w') as f:
    for res in res_python:
        f.write("%s " %res)
