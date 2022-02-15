import numpy as np
from Cell_funcs import *
from ENERGY.adhesions import *
from ENERGY.interaction import *
import numba as nb
import math

@nb.jit(nopython = True, nogil = True)
def tension_overlap_energy(cells, cparams, num_i, Ovlaps,
                E_const, h):

    oval_energy = 0.0

    x1, y1 = cells[3*num_i],cells[3*num_i+1]
    a1, b1 = cparams[4*num_i],cparams[4*num_i+1]
    theta_i = cparams[4*num_i+2] + cells[3*num_i+2]

    if np.any(Ovlaps[num_i] == 1):

        # Ellipse
        ells_i = create_ellipse((x1, y1), (a1, b1), theta_i)

        ind = np.arange(cells.shape[0]//3)[Ovlaps[num_i] == 1]

        for num_j in ind:
            x2, y2 = cells[3*num_j], cells[3*num_j+1]
            a2, b2 = cparams[4*num_j], cparams[4*num_j+1]
            theta_j = cparams[4*num_j+2] + cells[3*num_j+2]

            # If Contraction -> decrease 'a'
            if cparams[4*num_j+3] > 0:
                a2 -= h
            else:
                # protrusion -> increase 'h'
                a2 += 100*h

            # Ellipse
            ells_j = create_ellipse((x2, y2), (a2, b2), theta_j)

            # Constants
            theta_c = np.arctan2(y2 - y1, x2 - x1)
            k_const = (2.+0.5*np.cos(2.*(theta_c-theta_i)))*(2.+0.5*np.cos(2.*(theta_c-theta_j)))

            #find the overlap area
            energy = pair_overlap_energy(ells_i, ells_j, 
                            k_const, E_const)

            oval_energy += energy

    return oval_energy

# Virtual contraction
@nb.jit(nopython = True, nogil = True)
def virtual_contraction(cparam_, Adh_, h):

    # Make copies of array
    Adh = Adh_.copy()
    cparam = cparam_.copy()

    # Update the contractions
    c = h/cparam[0]

    # Update a
    cparam[0] -= h

    ind = np.arange(Adh.shape[0])[np.logical_and(
        Adh[:, 0] != -1e8, Adh[:, 1] != -1e8)]
    x, y = Adh[ind, 0], Adh[ind, 1]
    theta = cparam[2]

    Adh[ind, 0] = (x - c*cos(theta)*cos(theta)*x -
                      c*sin(theta)*cos(theta)*y)
    Adh[ind, 1] = (y - c*sin(theta)*cos(theta)*x -
                      c*sin(theta)*sin(theta)*y)

    return cparam, Adh

# Virtual extension
@nb.jit(nopython = True, nogil = True)
def virtual_extension(cparam_, Adh_, h):

    # Make copies of array
    Adh = Adh_.copy()
    cparam = cparam_.copy()

    # Update the extension
    c = h/cparam[0]

    # Update a
    cparam[0] += h

    ind = np.arange(Adh.shape[0])[np.logical_and(
        Adh[:, 0] != -1e8, Adh[:, 1] != -1e8)]
    x, y = Adh[ind, 0], Adh[ind, 1]
    theta = cparam[2]

    Adh[ind, 0] = (x + c*cos(theta)*cos(theta)*x +
                      c*sin(theta)*cos(theta)*y)
    Adh[ind, 1] = (y + c*sin(theta)*cos(theta)*x +
                      c*sin(theta)*sin(theta)*y)

    return cparam, Adh

#Compute tension
@nb.jit(nopython = True, nogil = True)
def compute_tension(cells_, num, cparams, Ovlaps_, Adh_, Adh0, E_const, k_s):

    cp_ = cparams.copy()
    cp = cparams.copy()
    cells = cells_.copy()
    Ovlaps = Ovlaps_.copy()
    Adh = Adh_.copy()

    step = 1e-3

    #Contraction - a - h decrease
    cp[4*num:4*num+4], Adh[num] = virtual_contraction(cp_[4*num:4*num+4],
                                Adh_[num], step)
    ce1 = tension_overlap_energy(cells, cp, num, Ovlaps,
                        E_const, step)
    ae1 = adhesion_energy(cells, num, Adh, Adh0, k_s)
    e1 = ce1 + ae1

    #Extension - a + h increase
    cp[4*num:4*num+4], Adh[num] = virtual_extension(cp_[4*num:4*num+4],
                                Adh_[num], step)
    ce2 = tension_overlap_energy(cells, cp, num, Ovlaps,
                        E_const, step)
    ae2 = adhesion_energy(cells, num, Adh, Adh0, k_s)
    e2 = ce2 + ae2

    #five point difference formula
    #(Force = -dU/dx)
    return -(e2 - e1)/(2*step)    