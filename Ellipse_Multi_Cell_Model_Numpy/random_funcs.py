import numpy as np
import numba as nb

from numpy.random import PCG64
from timeit import timeit

bit_gen = PCG64()
next_d = bit_gen.cffi.next_double
state = bit_gen.cffi.state_address

#Random double values
@nb.jit(nopython = True, nogil = True)
def uniform_double(low, high, n):

    if (high <= low):
        return None

    arr = np.empty(n)
    for i in range(n):
        arr[i] = low + (high - low)*next_d(state)

    return arr

@nb.jit(nopython = True, nogil = True)
def normals(n):
    out = np.empty(n)
    for i in range((n + 1) // 2):
        x1 = 2.0 * next_d(state) - 1.0
        x2 = 2.0 * next_d(state) - 1.0
        r2 = x1 * x1 + x2 * x2
        while r2 >= 1.0 or r2 == 0.0:
            x1 = 2.0 * next_d(state) - 1.0
            x2 = 2.0 * next_d(state) - 1.0
            r2 = x1 * x1 + x2 * x2
        f = np.sqrt(-2.0 * np.log(r2) / r2)
        out[2 * i] = f * x1
        if 2 * i + 1 < n:
            out[2 * i + 1] = f * x2
    return out

@nb.jit(nopython = True, nogil = True)
def gaussian(mu, sigma, n):

    out = normals(n)
    return (sigma*out + mu)