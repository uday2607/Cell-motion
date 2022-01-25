import numpy as np
from Cell_funcs import *
from Adh_funcs import *
import numba as nb

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

@nb.jit(nopython = True, nogil = True)
def pair_overlap_energy(ell_i_out, ell_i_in, ell_j_out, ell_j_in, k_const,
                    k_out_out, k_in_out, k_in_in):

    oval_energy = 0.0

    #find the overlap area
    area = ell_ell_area(ell_i_out, ell_j_out)
    if area > 0.0:
        oval_energy -= 0.5*k_out_out*k_const*area*area
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
def overlap_energy(cells, cparams, num_i, Ovlaps,
                k_out_out, k_in_out, k_in_in):

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

            # Ellipse
            ell_j_out, ell_j_in = create_ellipse((x2, y2), (a2, b2), theta_j)

            # Constants
            theta_c = np.arctan2(y2 - y1, x2 - x1)
            k_const = (1.+0.5*np.cos(2.*(theta_c-theta_i)))*(1.+0.5*np.cos(2.*(theta_c-theta_j)))

            #find the overlap area
            energy = pair_overlap_energy(ell_i_out, ell_i_in, ell_j_out, ell_j_in,
                        k_const, k_out_out, k_in_out, k_in_in)

            oval_energy += energy

    return oval_energy

@nb.jit(nopython = True, nogil = True)
def total_overlap_energy(cells, cparams, Ovlaps,
                k_out_out, k_in_out, k_in_in):

    total_energy = 0

    #Find the overlaps
    Ovlaps = find_overlaps(cells, cparams, Ovlaps)
    Eovlaps = np.zeros(Ovlaps.shape)

    for num_i in range(cells.shape[0]//3):

        #Cell coordinates
        x1, y1 = cells[3*num_i],cells[3*num_i+1]
        a1, b1 = cparams[4*num_i],cparams[4*num_i+1]
        theta_i = cparams[4*num_i+2] + cells[3*num_i+2]

        if np.any(Ovlaps[num_i] == 1):
            ind = np.arange(cells.shape[0]//3)[Ovlaps[num_i] == 1]

            # Ellipse
            ell_i_out, ell_i_in = create_ellipse((x1, y1), (a1, b1),
                                        theta_i)

            for num_j in ind:
                x2, y2 = cells[3*num_j],cells[3*num_j+1]
                a2, b2 = cparams[4*num_j],cparams[4*num_j+1]
                theta_j = cparams[4*num_j+2] + cells[3*num_j+2]

                # Ellipse
                ell_j_out, ell_j_in = create_ellipse((x2, y2), (a2, b2),
                                        theta_j)

                # Constants
                theta_c = np.arctan2(y2 - y1, x2 - x1)
                k_const = (1.+0.5*np.cos(2.*(theta_c-theta_i)))*(1.+0.5*np.cos(2.*(theta_c-theta_j)))

                if Eovlaps[num_i, num_j] == 0:
                    E = pair_overlap_energy(ell_i_out, ell_i_in, ell_j_out, ell_j_in,
                            k_const, k_out_out, k_in_out, k_in_in)
                    Eovlaps[num_i, num_j] = E
                    Eovlaps[num_i, num_j] = E

    total_energy += Eovlaps.sum()

    return total_energy

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
def compute_tension(cells_, num, cparams, Ovlaps_, Adh_, Adh0, k_out_out,
                    k_in_out, k_in_in, k_s):

    cp_ = cparams.copy()
    cp = cparams.copy()
    cells = cells_.copy()
    Ovlaps = Ovlaps_.copy()
    Adh = Adh_.copy()

    step = 1e-3

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
    return -(e1 - 8*e2 + 8*e3 - e4)/(12*step)

#Total Energy
@nb.jit(nopython = True, nogil = True)
def total_energy(cells, cparams, Ovlaps, Adh, Adh0,
                k_s, k_out_out, k_in_out, k_in_in, Y_MAX, Y_MIN):

    total_energy = 0

    #Find the overlaps
    Ovlaps = find_overlaps(cells, cparams, Ovlaps)
    Eovlaps = np.zeros(Ovlaps.shape)

    for num_i in range(cells.shape[0]//3):

        #Cell coordinates
        x1, y1 = cells[3*num_i],cells[3*num_i+1]
        a1, b1 = cparams[4*num_i],cparams[4*num_i+1]
        theta_i = cparams[4*num_i+2] + cells[3*num_i+2]

        # Ellipse
        ell_i_out, ell_i_in = create_ellipse((x1, y1), (a1, b1),
                                    theta_i)

        # constrained energy
        total_energy += 1e6*np.any(ell_i_out.T[1] > Y_MAX)
        total_energy += 1e6*np.any(ell_i_out.T[1] < Y_MIN)

        #Adhesion energy
        total_energy += adhesion_energy(cells, num_i, Adh, Adh0, k_s)

        if np.any(Ovlaps[num_i] == 1):
            ind = np.arange(cells.shape[0]//3)[Ovlaps[num_i] == 1]

            for num_j in ind:
                x2, y2 = cells[3*num_j],cells[3*num_j+1]
                a2, b2 = cparams[4*num_j],cparams[4*num_j+1]
                theta_j = cparams[4*num_j+2] + cells[3*num_j+2]

                # Ellipse
                ell_j_out, ell_j_in = create_ellipse((x2, y2), (a2, b2),
                                        theta_j)

                # Constants
                theta_c = np.arctan2(y2 - y1, x2 - x1)
                k_const = (1.+0.5*np.cos(2.*(theta_c-theta_i)))*(1.+0.5*np.cos(2.*(theta_c-theta_j)))

                if Eovlaps[num_i, num_j] == 0:
                    E = pair_overlap_energy(ell_i_out, ell_i_in, ell_j_out, ell_j_in,
                            k_const, k_out_out, k_in_out, k_in_in)
                    Eovlaps[num_i, num_j] = E
                    Eovlaps[num_j, num_i] = E

    total_energy += Eovlaps.sum()

    return total_energy

@nb.jit(nopython = True, nogil = True)
def total_energy_gradient(cells, cparams, Ovlaps, Adh, Adh0,
                k_s, k_out_out, k_in_out, k_in_in, Y_MAX, Y_MIN):

    #Vary every parameter by some small step and find
    #the partial differential
    grad = np.empty(cells.shape[0])
    e2 = 0.0
    e1 = 0.0
    step_x = 0.005
    step_th = 0.001
    ell_in = np.zeros((cells.shape[0]//3, 45, 2))
    ell_out = np.zeros((cells.shape[0]//3, 45, 2))
    adh_energy = np.zeros(cells.shape[0]//3)
    overlap_energy = np.zeros((cells.shape[0]//3, cells.shape[0]//3))

    #Ellipses
    for num_i in range(cells.shape[0]//3):
        #Cell coordinates
        x1, y1 = cells[3*num_i],cells[3*num_i+1]
        a1, b1 = cparams[4*num_i],cparams[4*num_i+1]
        theta_i = cparams[4*num_i+2] + cells[3*num_i+2]

        ell_i_out, ell_i_in = create_ellipse((x1, y1), (a1, b1), theta_i)
        ell_out[num_i] = ell_i_out
        ell_in[num_i] = ell_i_in

    # Compute the energy once to populate all the arrays
    # Find the overlaps
    Ovlaps = find_overlaps(cells, cparams, Ovlaps)

    for num_i in range(cells.shape[0]//3):

        #Adhesion energy
        adh_energy[num_i] = adhesion_energy(cells, num_i, Adh, Adh0, k_s)

        #Cell 1 coordinates
        x1, y1 = cells[3*num_i],cells[3*num_i+1]
        theta_i = cparams[4*num_i+2] + cells[3*num_i+2]

        if np.any(Ovlaps[num_i] == 1):
            ind = np.arange(cells.shape[0]//3)[Ovlaps[num_i] == 1]

            for num_j in ind:
                #Cell 2 cooridnates
                x2, y2 = cells[3*num_j],cells[3*num_j+1]
                theta_j = cparams[4*num_j+2] + cells[3*num_j+2]

                # Constants
                theta_c = np.arctan2(y2 - y1, x2 - x1)
                k_const = (1.+0.5*np.cos(2.*(theta_c-theta_i)))*(1.+0.5*np.cos(2.*(theta_c-theta_j)))

                if overlap_energy[num_i, num_j] == 0:
                    E = pair_overlap_energy(ell_out[num_i], ell_in[num_i],
                            ell_out[num_j], ell_in[num_j],
                            k_const, k_out_out, k_in_out, k_in_in)
                    overlap_energy[num_i, num_j] = E
                    overlap_energy[num_j, num_i] = E

    # Now compute the gradient
    for i in range(cells.shape[0]):

        # Make temporary copies
        temp_cell_i = cells[i]
        temp_ad_energy = adh_energy[i//3]
        temp_overlap_energy = overlap_energy.copy()
        temp_ovlap_ind = Ovlaps.copy()
        temp_ell_in = ell_in[i//3]
        temp_ell_out = ell_out[i//3]

        # change the variable
        if i % 3 != 2:
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

        ell_i_out, ell_i_in = create_ellipse((x1, y1), (a1, b1), theta_i)
        ell_out[num_i] = ell_i_out
        ell_in[num_i] = ell_i_in

        # constrained energy
        e2 += 1e6*np.any(ell_out[num_i].T[1] > Y_MAX)
        e2 += 1e6*np.any(ell_in[num_i].T[1] < Y_MIN)

        #Adhesion energy
        adh_energy[num_i] = adhesion_energy(cells, num_i, Adh, Adh0, k_s)

        #Overlap energy
        if np.any(Ovlaps[num_i] == 1):
            ind = np.arange(cells.shape[0]//3)[Ovlaps[num_i] == 1]

            for num_j in ind:
                #Cell 2 cooridnates
                x2, y2 = cells[3*num_j],cells[3*num_j+1]
                theta_j = cparams[4*num_j+2] + cells[3*num_j+2]

                # Constants
                theta_c = np.arctan2(y2 - y1, x2 - x1)
                k_const = (1.+0.5*np.cos(2.*(theta_c-theta_i)))*(1.+0.5*np.cos(2.*(theta_c-theta_j)))

                E = pair_overlap_energy(ell_out[num_i], ell_in[num_i],
                        ell_out[num_j], ell_in[num_j],
                        k_const, k_out_out, k_in_out, k_in_in)
                overlap_energy[num_i, num_j] = E
                overlap_energy[num_j, num_i] = E

        e2 = sum(adh_energy)
        e2 += overlap_energy.sum()

        #revert to original
        cells[i] = temp_cell_i
        adh_energy[i//3] = temp_ad_energy
        overlap_energy = temp_overlap_energy.copy()
        Ovlaps = temp_ovlap_ind.copy()
        ell_in[i//3] = temp_ell_in
        ell_out[i//3] = temp_ell_out

        # change the variable
        if i % 3 != 2:
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

        ell_i_out, ell_i_in = create_ellipse((x1, y1), (a1, b1), theta_i)
        ell_out[num_i] = ell_i_out
        ell_in[num_i] = ell_i_in

        # constrained energy
        e1 += 1e6*np.any(ell_out[num_i].T[1] > Y_MAX)
        e1 += 1e6*np.any(ell_in[num_i].T[1] < Y_MIN)

        #Adhesion energy
        adh_energy[num_i] = adhesion_energy(cells, num_i, Adh, Adh0, k_s)

        #Overlap energy
        if np.any(Ovlaps[num_i] == 1):
            ind = np.arange(cells.shape[0]//3)[Ovlaps[num_i] == 1]

            for num_j in ind:
                #Cell 2 cooridnates
                x2, y2 = cells[3*num_j],cells[3*num_j+1]
                theta_j = cparams[4*num_j+2] + cells[3*num_j+2]

                # Constants
                theta_c = np.arctan2(y2 - y1, x2 - x1)
                k_const = (1.+0.5*np.cos(2.*(theta_c-theta_i)))*(1.+0.5*np.cos(2.*(theta_c-theta_j)))

                E = pair_overlap_energy(ell_out[num_i], ell_in[num_i],
                        ell_out[num_j], ell_in[num_j],
                        k_const, k_out_out, k_in_out, k_in_in)
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
        ell_in[i//3] = temp_ell_in
        ell_out[i//3] = temp_ell_out

        # change the variable
        if i % 3 != 2:
            grad[i] = 0.5*(e2 - e1)/step_x
        else:
            grad[i] = 0.5*(e2 - e1)/step_th

    return grad
