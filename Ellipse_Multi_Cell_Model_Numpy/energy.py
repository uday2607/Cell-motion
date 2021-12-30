import numpy as np
from Cell_funcs import *
from Adh_funcs import *
import numba as nb

@nb.jit(nopython = True, nogil = True)
def adhesion_energy(cells, num, Adh, Adh0, k_s):

    ad_energy = 0
    fx = 0
    fy = 0

    ind = np.arange(Adh.shape[1])[np.logical_and(
        Adh[num, :, 0] != -1e8, Adh[num, :, 1] != -1e8)]

    dtheta = cells[3*num+2]

    for i in ind:
        fx = -(np.cos(dtheta)*Adh[num, i, 0]-np.sin(dtheta)*Adh[num, i, 1] +
                cells[3*num] - Adh0[num, i, 0])
        fy = -(np.sin(dtheta)*Adh[num, i, 0]+np.cos(dtheta)*Adh[num, i, 1] +
                cells[3*num+1] - Adh0[num, i, 1])

        ad_energy += 0.5*k_s*(fx**2 + fy**2)

    #ad_energy = 0.5*k_s*np.sum(fx**2 + fy**2)

    return ad_energy

@nb.jit(nopython = True, nogil = True)
def total_adhesion_energy(cells, Adh, Adh0, k_s):

    ad_energy = 0

    for i in range(cells.shape[0]//3):

        ad_energy += adhesion_energy(cells, i, Adh, Adh0, k_s)

    return ad_energy

@nb.jit(nopython = True, nogil = True)
def overlap_energy(cells, cparams, num, Ovlaps,
                    k_out_out, k_in_out, k_in_in):

    oval_energy = 0

    #cells overlap
    if np.any(Ovlaps[num] != -1e8):

        ell_i_out = create_ellipse((cells[3*num],cells[3*num+1]), (
                        cparams[4*num],cparams[4*num+1]),
                        cparams[4*num+2]+cells[3*num+2])
        ell_i_in = create_ellipse((cells[3*num],cells[3*num+1]), (0.5*
                        cparams[4*num],0.5*cparams[4*num+1]),
                        cparams[4*num+2]+cells[3*num+2])

        ind = np.arange(Ovlaps.shape[0])[Ovlaps[num] == 1]

        #for each cell, find the overlap energy
        for j in ind:
            ell_j_out = create_ellipse((cells[3*j],cells[3*j+1]), (
                            cparams[4*j],cparams[4*j+1]),
                            cparams[4*j+2]+cells[3*j+2])
            ell_j_in = create_ellipse((cells[3*j],cells[3*j+1]), (0.5*
                            cparams[4*j],0.5*cparams[4*j+1]),
                            cparams[4*j+2]+cells[3*j+2])

            #find the overlap area
            intersection = ellipse_intersection(ell_i_out, ell_j_out)
            if intersection.shape[0] != 0:
                oval_energy -= 0.5*k_out_out*polygon_area(intersection)
            intersection = ellipse_intersection(ell_i_in, ell_j_out)
            if intersection.shape[0] != 0:
                oval_energy += 0.5*k_in_out*polygon_area(intersection)
            intersection = ellipse_intersection(ell_i_out, ell_j_in)
            if intersection.shape[0] != 0:
                oval_energy += 0.5*k_in_out*polygon_area(intersection)
            intersection = ellipse_intersection(ell_i_in, ell_j_in)
            if intersection.shape[0] != 0:
                oval_energy -= 0.5*k_in_in*polygon_area(intersection)

    return oval_energy

@nb.jit(nopython = True, nogil = True)
def total_oval_energy(cells, cparams, Ovlaps,
                    k_out_out, k_in_out, k_in_in):

    total_oval_energy = 0

    #Find the overlaps
    Ovlaps = find_overlaps(cells, cparams, Ovlaps)
    Eovlaps = np.zeros(Ovlaps.shape)

    for num in range(cells.shape[0]//3):
        if np.any(Ovlaps[num]):
            ind = np.arange(cells.shape[0]//3)[Ovlaps[num] == 1]
            for i in ind:
                if Eovlaps[num, i] == 0:
                    E = overlap_energy(cells, cparams, num, Ovlaps,
                                k_out_out, k_in_out, k_in_in)
                    Eovlaps[num, i] = E
                    Eovlaps[i, num] = E
            total_oval_energy += np.sum(Eovlaps[num])

    return total_oval_energy

# Virtual contraction
@nb.jit(nopython = True, nogil = True)
def virtual_contraction(cparam_, Adh_, h):

    # Make copies of array
    Adh = np.array(Adh_)
    cparam = np.array(cparam_)

    # Update the contractions
    c = h/cparam[0]

    # Update a
    cparam[0] -= h

    ind = np.arange(Adh.shape[0])[np.logical_and(
        Adh[:, 0] != -1e8, Adh[:, 1] != -1e8)]
    x, y = Adh[ind, 0], Adh[ind, 1]
    theta = (cparam[2])

    Adh[ind, 0] = (x - c*cos(theta)*cos(theta)*x -
                      c*sin(theta)*cos(theta)*y)
    Adh[ind, 1] = (y - c*sin(theta)*cos(theta)*x -
                      c*sin(theta)*sin(theta)*y)

    return cparam, Adh

# Virtual extension
@nb.jit(nopython = True, nogil = True)
def virtual_extension(cparam_, Adh_, h):

    # Make copies of array
    Adh = np.array(Adh_)
    cparam = np.array(cparam_)

    # Update the extension
    c = h/cparam[0]

    # Update a
    cparam[0] += h

    ind = np.arange(Adh.shape[0])[np.logical_and(
        Adh[:, 0] != -1e8, Adh[:, 1] != -1e8)]
    x, y = Adh[ind, 0], Adh[ind, 1]
    theta = (cparam[2])

    Adh[ind, 0] = (x + c*cos(theta)*cos(theta)*x +
                      c*sin(theta)*cos(theta)*y)
    Adh[ind, 1] = (y + c*sin(theta)*cos(theta)*x +
                      c*sin(theta)*sin(theta)*y)

    return cparam, Adh

#Compute tension
@nb.jit(nopython = True, nogil = True)
def compute_tension(cells_, num, cparams, Ovlaps_, Adh_, Adh0, k_out_out,
                    k_in_out, k_in_in, k_s):

    cp_ = np.array(cparams)
    cp = np.array(cparams)
    cells = np.array(cells_)
    Ovlaps = np.array(Ovlaps_)
    Adh = np.array(Adh_)

    step = 1e-4

    #Extension - a + h increase
    cp[4*num:4*num+4], Adh[num] = virtual_extension(cp_[4*num:4*num+4],
                                Adh_[num], step)
    ce3 = overlap_energy(cells, cp, num, Ovlaps,
                        k_out_out, k_in_out, k_in_in)
    ae3 = adhesion_energy(cells, num, Adh, Adh0, k_s)
    e3 = ce3 + ae3

    #Extension - a + 2h increase
    cp[4*num:4*num+4], Adh[num] = virtual_extension(cp[4*num:4*num+4],
                                Adh[num], step)
    ce4 = overlap_energy(cells, cp, num, Ovlaps,
                        k_out_out, k_in_out, k_in_in)
    ae4 = adhesion_energy(cells, num, Adh, Adh0, k_s)
    e4 = ce4 + ae4

    #Contraction - a - h decrease
    cp[4*num:4*num+4], Adh[num] = virtual_contraction(cp_[4*num:4*num+4],
                                Adh_[num], step)
    ce2 = overlap_energy(cells, cp, num, Ovlaps,
                        k_out_out, k_in_out, k_in_in)
    ae2 = adhesion_energy(cells, num, Adh, Adh0, k_s)
    e2 = ce2 + ae2

    #Contraction - a - 2h decrease
    cp[4*num:4*num+4], Adh[num] = virtual_contraction(cp[4*num:4*num+4],
                                Adh[num], step)
    ce1 = overlap_energy(cells, cp, num, Ovlaps,
                        k_out_out, k_in_out, k_in_in)
    ae1 = adhesion_energy(cells, num, Adh, Adh0, k_s)
    e1 = ce1 + ae1

    #five point difference formula
    #(Force = -dU/dx)
    return -1*(e1 - 8*e2 + 8*e3 - e4)/(12*step)

#Total Energy
@nb.jit(nopython = True, nogil = True)
def total_energy(cells, cparams, Ovlaps, Adh, Adh0,
                k_s, k_out_out, k_in_out, k_in_in):

    return (total_adhesion_energy(cells, Adh, Adh0, k_s) +
            total_oval_energy(cells, cparams, Ovlaps,
                                k_out_out, k_in_out, k_in_in))

@nb.jit(nopython = True, nogil = True)
def total_energy_gradient(cells, cparams, Ovlaps, Adh, Adh0,
                k_s, k_out_out, k_in_out, k_in_in):

    #Vary every parameter by some small step and find
    #the partial differential
    grad = np.empty(cells.shape)
    cells0 = np.array(cells)
    step = 1e-4

    for i in range(cells.shape[0]):
        cells0[i] = cells[i] - 2*step
        e1 = total_energy(cells0, cparams, Ovlaps, Adh, Adh0,
                k_s, k_out_out, k_in_out, k_in_in)

        cells0[i] = cells[i] - step
        e2 = total_energy(cells0, cparams, Ovlaps, Adh, Adh0,
                k_s, k_out_out, k_in_out, k_in_in)

        cells0[i] = cells[i] + step
        e3 = total_energy(cells0, cparams, Ovlaps, Adh, Adh0,
                k_s, k_out_out, k_in_out, k_in_in)

        cells0[i] = cells[i] + 2*step
        e4 = total_energy(cells0, cparams, Ovlaps, Adh, Adh0,
                k_s, k_out_out, k_in_out, k_in_in)

        grad[i] = (e1 - 8*e2 + 8*e3 - e4)/(12*step)

    return grad

@nb.jit(nopython = True, nogil = True)
def total_energy_and_grad(cells, cparams, Ovlaps, Adh, Adh0,
                k_s, k_out_out, k_in_out, k_in_in):

    return np.array([total_energy(cells, cparams, Ovlaps, Adh, Adh0,
                k_s, k_out_out, k_in_out, k_in_in),
                total_energy_and_grad(cells, cparams, Ovlaps, Adh, Adh0,
                k_s, k_out_out, k_in_out, k_in_in)])
