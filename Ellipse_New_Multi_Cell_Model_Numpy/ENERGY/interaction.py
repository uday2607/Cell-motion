import numpy as np
from Cell_funcs import *
import numba as nb
import math

@nb.jit(nopython = True, nogil = True)
def pair_overlap_energy(ells_i, ells_j, k_const,
                    E_const):

    oval_energy = 0.0

    #find the overlap area
    #out - out ellipses
    area = ell_ell_area(ells_i[2], ells_j[2])
    if area > 0.0:
        oval_energy += 0.5*E_const[2, 2]*k_const*area*area
        
        #mid - mid ellipses
        area = ell_ell_area(ells_i[1], ells_j[1])
        if area > 0.0:
            oval_energy += 0.5*E_const[1, 1]*k_const*area*area

            #mid - in ellipses
            area = ell_ell_area(ells_i[1], ells_j[0])
            if area > 0.0:
                oval_energy += 0.5*E_const[1, 0]*k_const*area*area

            #in - mid ellipses
            area = ell_ell_area(ells_i[0], ells_j[1])
            if area > 0.0:
                oval_energy += 0.5*E_const[0, 1]*k_const*area*area  

                #in - in ellipses
                area = ell_ell_area(ells_i[0], ells_j[0])
                if area > 0.0:
                    oval_energy += 0.5*E_const[0, 0]*k_const*area*area  

    return oval_energy

@nb.jit(nopython = True, nogil = True)
def total_overlap_energy(cells, cparams, Ovlaps,
                E_const):

    total_energy = 0

    #Find the overlaps
    Ovlaps = find_overlaps(cells, cparams, Ovlaps)
    Eovlaps = np.zeros(Ovlaps.shape)

    for num_i in range(cells.shape[0]//3):

        #Cell coordinates
        x1, y1 = cells[3*num_i], cells[3*num_i+1]
        a1, b1 = cparams[4*num_i], cparams[4*num_i+1]
        theta_i = cparams[4*num_i+2] + cells[3*num_i+2]
        theta_i -= 2.0*pi*math.floor((theta_i + pi)*(1.0/(2.0*pi)))

        if np.any(Ovlaps[num_i] == 1):
            ind = np.arange(cells.shape[0]//3)[Ovlaps[num_i] == 1]

            # Ellipse
            ells_i = create_ellipse((x1, y1), (a1, b1),
                                        theta_i)

            for num_j in ind:
                x2, y2 = cells[3*num_j],cells[3*num_j+1]
                a2, b2 = cparams[4*num_j],cparams[4*num_j+1]
                theta_j = cparams[4*num_j+2] + cells[3*num_j+2]
                theta_j -= 2.0*pi*math.floor((theta_j + pi)*(1.0/(2.0*pi)))

                # Ellipse
                ells_j = create_ellipse((x2, y2), (a2, b2),
                                        theta_j)

                # Constants
                theta_c = np.arctan2(y2 - y1, x2 - x1)
                k_const = (2.+0.5*np.cos(2.*(theta_c-theta_i)))*(2.+0.5*np.cos(2.*(theta_c-theta_j)))

                if Eovlaps[num_i, num_j] == 0:
                    E = pair_overlap_energy(ells_i, ells_j,
                            k_const, E_const)
                    Eovlaps[num_i, num_j] = E
                    Eovlaps[num_i, num_j] = E

    total_energy += Eovlaps.sum()

    return total_energy