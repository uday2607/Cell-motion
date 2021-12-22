import numpy as np
from cell import *
from adhesions import *
import numba as nb

def adhesion_energy(cells, Adh, Adh0, k_s):

    ad_energy = 0

    for i in range(cells.shape[0]//3):
        ind = np.arange(Adh.shape[1])[np.all(Adh[i]!=-1e8,axis=1)]
        dtheta = pi/180*cells[3*i+2]
        fx = -((cos(dtheta)*Adh[i, ind, 0]-sin(dtheta)*Adh[i, ind, 1] +
                cells[3*i] - Adh0[i, ind, 0]))
        fy = -((sin(dtheta)*Adh[i, ind, 0]+cos(dtheta)*Adh[i, ind, 1] +
                cells[3*i+1] - Adh0[i, ind, 1]))

        ad_energy += 0.5*k_s*np.sum(fx**2 + fy**2)

    return ad_energy

def overlap_energy(cells, cparams, num, Ovlaps,
                    k_out_out, k_in_out, k_in_in):

    oval_energy = 0

    #cells overlap
    if np.any(Ovlaps[num] != -1e8):

        ell_i_out = create_ellipse((cells[3*num],cells[3*num+1]), (
                        cparams[4*num],cparams[4*num+1]), cparams[4*num+2])
        ell_i_in = create_ellipse((cells[3*num],cells[3*num+1]), (0.5*
                        cparams[4*num],0.5*cparams[4*num+1]), cparams[4*num+2])

        ind = np.arange(Ovlaps.shape[0])[Ovlaps[num] == 1]

        #for each cell, find the overlap energy
        for j in ind:
            ell_j_out = create_ellipse((cells[3*j],cells[3*j+1]), (
                            cparams[4*j],cparams[4*j+1]), cparams[4*j+2])
            ell_j_in = create_ellipse((cells[3*j],cells[3*j+1]), (0.5*
                            cparams[4*j],0.5*cparams[4*j+1]), cparams[4*j+2])

            #find the overlap area
            oval_energy -= 0.5*k_out_out*(ell_i_out.intersection(ell_j_out).area)
            oval_energy += 0.5*k_in_out*(ell_i_out.intersection(ell_j_in).area)
            oval_energy += 0.5*k_in_out*(ell_i_in.intersection(ell_j_out).area)
            oval_energy += 0.5*k_in_in*(ell_i_in.intersection(ell_j_in).area)

    return oval_energy

def total_oval_energy(cells, cparams, Ovlaps,
                    k_out_out, k_in_out, k_in_in):

    total_oval_energy = 0

    for num in range(cells.shape[0]//3):
        total_oval_energy += overlap_energy(cells, cparams, num, Ovlaps,
                            k_out_out, k_in_out, k_in_in)

    return total_oval_energy

def compute_tension(cells, num, cparams, Ovlaps, k_out_out,
                    k_in_out, k_in_in):

    cp = cparams.copy()

    #change a
    cp[4*num] = 1.5*cparams[4*num]
    e1 = overlap_energy(cells, cp, num, Ovlaps,
                        k_out_out, k_in_out, k_in_in)

    cp[4*num] = 0.5*cparams[4*num]
    e2 = overlap_energy(cells, cp, num, Ovlaps,
                        k_out_out, k_in_out, k_in_in)

    return (e2 - e1)/cparams[4*num]

def total_energy(cells, cparams, Ovlaps, Adh, Adh0,
                k_s, k_out_out, k_in_out, k_in_in):

    return (adhesion_energy(cells, Adh, Adh0, k_s) +
            total_oval_energy(cells, cparams, Ovlaps,
                                k_out_out, k_in_out, k_in_in))
