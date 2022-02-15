import numpy as np
from Cell_funcs import *
from ENERGY.adhesions import *
import numba as nb

# Of inner ellipses only
@nb.jit(nopython = True, nogil = True)
def pair_overlap_energy(ell_i_out, ell_i_in, ell_j_out, ell_j_in, k_const,
                    k_out_out, k_in_out, k_in_in):

    oval_energy = 0.0

    #find the overlap area
    area = ell_ell_area(ell_i_out, ell_j_out)
    if area > 0.0:
        oval_energy -= 0.0 # No force from out-out ellipses
        #in out
        area = ell_ell_area(ell_i_in, ell_j_out)
        if area > 0.0:
            oval_energy += 0.5*k_in_out*k_const*area*area
        #out in
        area = ell_ell_area(ell_i_out, ell_j_in)
        if area > 0.0:
            oval_energy += 0.5*k_in_out*k_const*area*area
            #in in
            area = ell_ell_area(ell_i_in, ell_j_in)
            if area > 0.0:
                oval_energy += 0.5*k_in_in*k_const*area*area

    return oval_energy

@nb.jit(nopython = True, nogil = True)
def force_overlap_energy(cells, cparams, num_i, Ovlaps,
                k_out_out, k_in_out, k_in_in, h):

    oval_energy = 0.0

    x1, y1 = cells[3*num_i],cells[3*num_i+1]
    a1, b1 = cparams[4*num_i],cparams[4*num_i+1]
    theta_i = cparams[4*num_i+2] + cells[3*num_i+2]

    if np.any(Ovlaps[num_i] == 1):

        # Ellipse
        ell_i_out, ell_i_in = create_ellipse((x1, y1), (a1, b1), theta_i)

        ind = np.arange(cells.shape[0]//3)[Ovlaps[num_i] == 1]

        for num_j in ind:
            x2, y2 = cells[3*num_j],cells[3*num_j+1]
            a2, b2 = cparams[4*num_j],cparams[4*num_j+1]
            theta_j = cparams[4*num_j+2] + cells[3*num_j+2]

            # If Contraction -> decrease 'a'
            if cparams[4*num_j+3] > 0:
                a2 -= h
            else:
                # protrusion -> increase 'h'
                a2 += 10*h
                #x2 += 10*h*cos(theta_j)
                #y2 += 10*h*sin(theta_j)

            # Ellipse
            ell_j_out, ell_j_in = create_ellipse((x2, y2), (a2, b2), theta_j)

            # Constants
            theta_c = np.arctan2(y2 - y1, x2 - x1)
            k_const = (2.+1.0*np.cos(2.*(theta_c-theta_i)))*(2.+1.0*np.cos(2.*(theta_c-theta_j)))

            #find the overlap area
            energy = pair_overlap_energy(ell_i_out, ell_i_in, ell_j_out, ell_j_in,
                        k_const, k_out_out, k_in_out, k_in_in)

            oval_energy += energy

    return oval_energy


# Virtual extension
@nb.jit(nopython = True, nogil = True)
def virtual_extension(cell_, cparam_, Adh_, h):

    # Make copies of array
    cell = cell_.copy()
    Adh = Adh_.copy()
    cparam = cparam_.copy()

    # Update the center of the cell
    theta = cparam[2]+cell[2]
    cell += 10*h*np.array([cos(theta), sin(theta)])

    # Update a
    cparam[0] += 10*h

    ind = np.arange(Adh.shape[0])[np.logical_and(
        Adh[:, 0] != -1e8, Adh[:, 1] != -1e8)]
    x, y = Adh[ind, 0], Adh[ind, 1]

    Adh[ind, 0] = x + cell_[0] - cell[0]
    Adh[ind, 1] = y + cell_[1] - cell[1]

    return cell, cparam, Adh

#Compute tension
@nb.jit(nopython = True, nogil = True)
def compute_force(cells_, num, cparams, Ovlaps_, Adh_, Adh0, k_out_out,
                    k_in_out, k_in_in, k_s):

    cp_ = cparams.copy()
    cp = cparams.copy()
    cells = cells_.copy()
    Ovlaps = Ovlaps_.copy()
    Adh = Adh_.copy()

    step = 1e-3

    ce1 = force_overlap_energy(cells, cp, num, Ovlaps,
                        k_out_out, k_in_out, k_in_in, step)
    ae1 = adhesion_energy(cells, num, Adh, Adh0, k_s)
    e1 = ce1 + ae1

    #Extension - a + h increase
    cells[3*num:3*num+3], cp[4*num:4*num+4], Adh[num] = virtual_extension(cells_[3*num:3*num+3],
                                cp_[4*num:4*num+4],
                                Adh_[num], step)
    ce2 = force_overlap_energy(cells, cp, num, Ovlaps,
                        k_out_out, k_in_out, k_in_in, step)
    ae2 = adhesion_energy(cells, num, Adh, Adh0, k_s)
    e2 = ce2 + ae2

    #Extension - a + 2h increase
    cells[3*num:3*num+3], cp[4*num:4*num+4], Adh[num] = virtual_extension(cells[3*num:3*num+3],
                                cp[4*num:4*num+4],
                                Adh[num], step)
    ce3 = force_overlap_energy(cells, cp, num, Ovlaps,
                        k_out_out, k_in_out, k_in_in, step)
    ae3 = adhesion_energy(cells, num, Adh, Adh0, k_s)
    e3 = ce3 + ae3

    #five point difference formula
    #(Force = -dU/dx)
    return (-3*e1 + 4*e2 - e3)/(2*step)