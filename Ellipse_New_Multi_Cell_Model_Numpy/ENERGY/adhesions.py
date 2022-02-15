import numpy as np
import numba as nb
import math

@nb.jit(nopython = True, nogil = True)
def adhesion_energy(cells, num, Adh, Adh0, k_s):

    ad_energy = 0.0
    fx = 0.0
    fy = 0.0
    cos_t = np.cos(cells[3*num+2])
    sin_t = np.sin(cells[3*num+2])

    ind = np.arange(Adh.shape[1])[np.logical_and(
        Adh[num, :, 0] != -1e8, Adh[num, :, 1] != -1e8)]

    for i in ind:
        fx = -(cos_t*Adh[num, i, 0]-sin_t*Adh[num, i, 1] +
                cells[3*num] - Adh0[num, i, 0])
        fy = -(sin_t*Adh[num, i, 0]+cos_t*Adh[num, i, 1] +
                cells[3*num+1] - Adh0[num, i, 1])

        ad_energy += (fx**2 + fy**2)

    return 0.5*k_s*ad_energy