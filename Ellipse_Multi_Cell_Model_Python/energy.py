import numpy as np
from cell import *
from adhesions import *
import numba as nb

def adhesion_energy(cells, cp, num, Adh, Adh0, k_s):

    ad_energy = 0

    ind = np.arange(Adh.shape[1])[np.all(Adh[num]!=-1e8,axis=1)]
    dtheta = pi/180*cells[3*num+2]
    fx = -((cos(dtheta)*Adh[num, ind, 0]-sin(dtheta)*Adh[num, ind, 1] +
            cells[3*num] - Adh0[num, ind, 0]))
    fy = -((sin(dtheta)*Adh[num, ind, 0]+cos(dtheta)*Adh[num, ind, 1] +
            cells[3*num+1] - Adh0[num, ind, 1]))

    ad_energy = 0.5*k_s*np.sum(fx**2 + fy**2)

    return ad_energy


def total_adhesion_energy(cells, cparams, Adh, Adh0, k_s):

    ad_energy = 0

    for i in range(cells.shape[0]//3):

        ad_energy += adhesion_energy(cells, cparams, i, Adh, Adh0, k_s)

    return ad_energy

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
            oval_energy -= 0.5*k_out_out*(ell_i_out.intersection(ell_j_out).area)
            if ell_i_out.boundary.intersects(ell_j_in.boundary):
                oval_energy += 0.5*k_in_out*(ell_i_out.intersection(ell_j_in).area)
            if ell_j_out.boundary.intersects(ell_i_in.boundary):
                oval_energy += 0.5*k_in_out*(ell_j_out.intersection(ell_i_in).area)
            if ell_i_in.boundary.intersects(ell_j_in.boundary):
                oval_energy += 0.5*k_in_in*(ell_i_in.intersection(ell_j_in).area)    

    return oval_energy, Ovlaps

def total_oval_energy(cells, cparams, Ovlaps,
                    k_out_out, k_in_out, k_in_in):

    total_oval_energy = 0

    #Find the overlaps
    Ovlaps = find_overlaps(cells, cparams, Ovlaps)

    for num in range(cells.shape[0]//3):
        total_oval_energy += overlap_energy(cells, cparams, num, Ovlaps,
                            k_out_out, k_in_out, k_in_in)[0]

    return total_oval_energy

# Virtual contraction
def virtual_contraction(cparam_, Adh_, h):

    # Make copies of array
    Adh = np.array(Adh_)
    cparam = np.array(cparam_)

    # Update the contractions
    c = h/cparam[0]

    # Update a
    cparam[0] -= h

    ind = np.arange(Adh.shape[0])[np.all(Adh!=-1e8,axis=1)]
    x, y = Adh[ind, 0], Adh[ind, 1]
    theta = pi/180*(cparam[2])

    Adh[ind, 0] = (x - c*cos(theta)*cos(theta)*x -
                      c*sin(theta)*cos(theta)*y)
    Adh[ind, 1] = (y - c*sin(theta)*cos(theta)*x -
                      c*sin(theta)*sin(theta)*y)

    return cparam, Adh

# Virtual extension
def virtual_extension(cparam_, Adh_, h):

    # Make copies of array
    Adh = np.array(Adh_)
    cparam = np.array(cparam_)

    # Update the extension
    c = h/cparam[0]

    # Update a
    cparam[0] += h

    ind = np.arange(Adh.shape[0])[np.all(Adh!=-1e8,axis=1)]
    x, y = Adh[ind, 0], Adh[ind, 1]
    theta = pi/180*(cparam[2])

    Adh[ind, 0] = (x + c*cos(theta)*cos(theta)*x +
                      c*sin(theta)*cos(theta)*y)
    Adh[ind, 1] = (y + c*sin(theta)*cos(theta)*x +
                      c*sin(theta)*sin(theta)*y)

    return cparam, Adh

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
    ce3, _ = overlap_energy(cells, cp, num, Ovlaps,
                        k_out_out, k_in_out, k_in_in)
    ae3 = adhesion_energy(cells, cp, num, Adh, Adh0, k_s)
    e3 = ce3 + ae3

    #Extension - a + 2h increase
    cp[4*num:4*num+4], Adh[num] = virtual_extension(cp[4*num:4*num+4], 
                                Adh[num], step)
    ce4, _ = overlap_energy(cells, cp, num, Ovlaps,
                        k_out_out, k_in_out, k_in_in)
    ae4 = adhesion_energy(cells, cp, num, Adh, Adh0, k_s)
    e4 = ce4 + ae4                 

    #Contraction - a - h decrease
    cp[4*num:4*num+4], Adh[num] = virtual_contraction(cp_[4*num:4*num+4], 
                                Adh_[num], step)
    ce2, _ = overlap_energy(cells, cp, num, Ovlaps,
                        k_out_out, k_in_out, k_in_in)
    ae2 = adhesion_energy(cells, cp, num, Adh, Adh0, k_s)
    e2 = ce2 + ae2

    #Contraction - a - 2h decrease
    cp[4*num:4*num+4], Adh[num] = virtual_contraction(cp[4*num:4*num+4], 
                                Adh[num], step)
    ce1, _ = overlap_energy(cells, cp, num, Ovlaps,
                        k_out_out, k_in_out, k_in_in)
    ae1 = adhesion_energy(cells, cp, num, Adh, Adh0, k_s)
    e1 = ce1 + ae1

    #five point difference formula
    return -1*(e1 - 8*e2 + 8*e3 - e4)/(12*step)

def total_energy(cells, cparams, Ovlaps, Adh, Adh0,
                k_s, k_out_out, k_in_out, k_in_in):

    #Find overlaps
    Ovlaps = find_overlaps(cells, cparams, Ovlaps)            

    return (total_adhesion_energy(cells, cparams, Adh, Adh0, k_s) +
            total_oval_energy(cells, cparams, Ovlaps,
                                k_out_out, k_in_out, k_in_in))
