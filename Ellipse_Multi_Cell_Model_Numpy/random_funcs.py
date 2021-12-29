import numpy as np
import numba as nb

from numpy.random import PCG64
from timeit import timeit

bit_gen = PCG64()
next_d = bit_gen.cffi.next_double
state_addr = bit_gen.cffi.state_address

#Random double values
@nb.jit(nopython = True, nogil = True)
def uniform_double(low, high, n):

    if (high <= low):
        return None

    arr = np.empty(n)
    for i in range(n):
        arr[i] = low + (high - low)*next_d(state_addr)

    return arr