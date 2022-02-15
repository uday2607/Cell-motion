import numpy as np
from Cell_funcs import *
from ENERGY.adhesions import *
from ENERGY.interaction import *
import numba as nb
import math

#Total Energy
@nb.jit(nopython = True, nogil = True, cache=True)
def total_energy(cells, cparams, Ovlaps, Adh, Adh0,
                k_s, E_const):

    total_energy = 0

    #Find the overlaps
    Ovlaps = find_overlaps(cells, cparams, Ovlaps)
    Eovlaps = np.zeros(Ovlaps.shape)

    for num_i in range(cells.shape[0]//3):

        #Adhesion energy
        total_energy += adhesion_energy(cells, num_i, Adh, Adh0, k_s)

        #Cell coordinates
        x1, y1 = cells[3*num_i],cells[3*num_i+1]
        a1, b1 = cparams[4*num_i],cparams[4*num_i+1]
        theta_i = cparams[4*num_i+2] + cells[3*num_i+2]
        #theta_i -= 2.0*pi*math.floor((theta_i + pi)*(1.0/(2.0*pi)))

        if np.any(Ovlaps[num_i] == 1):
            ind = np.arange(cells.shape[0]//3)[Ovlaps[num_i] == 1]

            # Ellipse
            ells_i = create_ellipse((x1, y1), (a1, b1),
                                        theta_i)

            for num_j in ind:
                x2, y2 = cells[3*num_j],cells[3*num_j+1]
                a2, b2 = cparams[4*num_j],cparams[4*num_j+1]
                theta_j = cparams[4*num_j+2] + cells[3*num_j+2]
                #theta_j -= 2.0*pi*math.floor((theta_j + pi)*(1.0/(2.0*pi)))

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
                    Eovlaps[num_j, num_i] = E

    total_energy += Eovlaps.sum()

    return total_energy

@nb.jit(nopython = True, nogil = True, cache=True)
def total_energy_gradient(cells, cparams, Ovlaps, Adh, Adh0,
                k_s, E_const):

    #Vary every parameter by some small step and find
    #the partial differential
    grad = np.empty(cells.shape[0])
    e2 = 0.0
    e1 = 0.0
    step_x = 1e-4
    step_th = 1e-4
    ells = np.zeros((cells.shape[0]//3, 3, 45, 2))
    adh_energy = np.zeros(cells.shape[0]//3)
    overlap_energy = np.zeros((cells.shape[0]//3, cells.shape[0]//3))

    #Ellipses
    for num_i in range(cells.shape[0]//3):
        #Cell coordinates
        x1, y1 = cells[3*num_i],cells[3*num_i+1]
        a1, b1 = cparams[4*num_i],cparams[4*num_i+1]
        theta_i = cparams[4*num_i+2] + cells[3*num_i+2]
        #theta_i -= 2.0*pi*math.floor((theta_i + pi)*(1.0/(2.0*pi)))

        ells_i = create_ellipse((x1, y1), (a1, b1), theta_i)
        ells[num_i] = ells_i

    # Compute the energy once to populate all the arrays
    # Find the overlaps
    Ovlaps = find_overlaps(cells, cparams, Ovlaps)

    for num_i in range(cells.shape[0]//3):

        #Adhesion energy
        adh_energy[num_i] = adhesion_energy(cells, num_i, Adh, Adh0, k_s)

        #Cell 1 coordinates
        x1, y1 = cells[3*num_i],cells[3*num_i+1]
        theta_i = cparams[4*num_i+2] + cells[3*num_i+2]
        #theta_i -= 2.0*pi*math.floor((theta_i + pi)*(1.0/(2.0*pi)))

        if np.any(Ovlaps[num_i] == 1):
            ind = np.arange(cells.shape[0]//3)[Ovlaps[num_i] == 1]

            for num_j in ind:
                #Cell 2 cooridnates
                x2, y2 = cells[3*num_j],cells[3*num_j+1]
                theta_j = cparams[4*num_j+2] + cells[3*num_j+2]
                #theta_j -= 2.0*pi*math.floor((theta_j + pi)*(1.0/(2.0*pi)))

                # Constants
                theta_c = np.arctan2(y2 - y1, x2 - x1)
                k_const = (2.+0.5*np.cos(2.*(theta_c-theta_i)))*(2.+0.5*np.cos(2.*(theta_c-theta_j)))

                if overlap_energy[num_i, num_j] == 0:
                    E = pair_overlap_energy(ells[num_i], ells[num_j],
                            k_const, E_const)
                    overlap_energy[num_i, num_j] = E
                    overlap_energy[num_j, num_i] = E

    # Now compute the gradient
    for i in range(cells.shape[0]):

        # Make temporary copies
        temp_cell_i = cells[i]
        temp_ad_energy = adh_energy[i//3]
        temp_overlap_energy = overlap_energy.copy()
        temp_ovlap_ind = Ovlaps.copy()
        temp_ells = ells[i//3].copy()

        # change the variable
        if i % 3 == 0:
            cells[i] += step_x
        elif i % 3 == 1:
            cells[i] += step_x
        else:
            cells[i] += step_th

        # Compute the ovlap indices
        Ovlaps = find_overlaps_ind(cells, cparams, Ovlaps, i//3)

        # Index of the cell
        #Cell coordinates
        num_i = i//3
        x1, y1 = cells[3*num_i],cells[3*num_i+1]
        a1, b1 = cparams[4*num_i],cparams[4*num_i+1]
        theta_i = cparams[4*num_i+2] + cells[3*num_i+2]
        #theta_i -= 2.0*pi*math.floor((theta_i + pi)*(1.0/(2.0*pi)))

        ells_i = create_ellipse((x1, y1), (a1, b1), theta_i)

        #Adhesion energy
        adh_energy[num_i] = adhesion_energy(cells, num_i, Adh, Adh0, k_s)

        #Overlap energy
        if np.any(Ovlaps[num_i] == 1):
            ind = np.arange(cells.shape[0]//3)[Ovlaps[num_i] == 1]

            for num_j in ind:
                #Cell 2 cooridnates
                x2, y2 = cells[3*num_j],cells[3*num_j+1]
                theta_j = cparams[4*num_j+2] + cells[3*num_j+2]
                #theta_j -= 2.0*pi*math.floor((theta_j + pi)*(1.0/(2.0*pi)))

                # Constants
                theta_c = np.arctan2(y2 - y1, x2 - x1)
                k_const = (2.+0.5*np.cos(2.*(theta_c-theta_i)))*(2.+0.5*np.cos(2.*(theta_c-theta_j)))

                E = pair_overlap_energy(ells_i, ells[num_j],
                        k_const, E_const)
                overlap_energy[num_i, num_j] = E
                overlap_energy[num_j, num_i] = E

        e2 = sum(adh_energy)
        e2 += overlap_energy.sum()

        #revert to original
        cells[i] = temp_cell_i
        adh_energy[i//3] = temp_ad_energy
        overlap_energy = temp_overlap_energy.copy()
        Ovlaps = temp_ovlap_ind.copy()
        ells[i//3] = temp_ells

        # change the variable
        if i % 3 == 0:
            cells[i] -= step_x
        elif i % 3 == 1:
            cells[i] -= step_x   
        else:
            cells[i] -= step_th

        # Compute the ovlap indices
        Ovlaps = find_overlaps_ind(cells, cparams, Ovlaps, i//3)

        # Index of the cell
        #Cell coordinates
        num_i = i//3
        x1, y1 = cells[3*num_i],cells[3*num_i+1]
        a1, b1 = cparams[4*num_i],cparams[4*num_i+1]
        theta_i = cparams[4*num_i+2] + cells[3*num_i+2]
        #theta_i -= 2.0*pi*math.floor((theta_i + pi)*(1.0/(2.0*pi)))

        ells_i = create_ellipse((x1, y1), (a1, b1), theta_i)

        #Adhesion energy
        adh_energy[num_i] = adhesion_energy(cells, num_i, Adh, Adh0, k_s)

        #Overlap energy
        if np.any(Ovlaps[num_i] == 1):
            ind = np.arange(cells.shape[0]//3)[Ovlaps[num_i] == 1]

            for num_j in ind:
                #Cell 2 cooridnates
                x2, y2 = cells[3*num_j],cells[3*num_j+1]
                theta_j = cparams[4*num_j+2] + cells[3*num_j+2]
                #theta_j -= 2.0*pi*math.floor((theta_j + pi)*(1.0/(2.0*pi)))

                # Constants
                theta_c = np.arctan2(y2 - y1, x2 - x1)
                k_const = (2.+0.5*np.cos(2.*(theta_c-theta_i)))*(2.+0.5*np.cos(2.*(theta_c-theta_j)))

                E = pair_overlap_energy(ells_i, ells[num_j],
                        k_const, E_const)
                overlap_energy[num_i, num_j] = E
                overlap_energy[num_j, num_i] = E

        e1 = sum(adh_energy)
        e1 += overlap_energy.sum()

        #revert to original
        # Make temporary copies
        cells[i] = temp_cell_i
        adh_energy[i//3] = temp_ad_energy
        overlap_energy = temp_overlap_energy.copy()
        Ovlaps = temp_ovlap_ind.copy()
        ells[i//3] = temp_ells

        # change the variable
        if i % 3 != 2:
            grad[i] = 0.5*(e2 - e1)/step_x
        else:
            grad[i] = 0.5*(e2 - e1)/step_th

    return grad/grad.shape[0]