import numpy as np

def adhesion_energy(cells, Adh, Adh0, k_s):

    ad_energy = 0

    for i in range(cells.shape[0]//3):
        ind = np.arange(Adh.shape[1])[Adh[i] != np.array([-1, -1])]
        theta = pi/180*cell[3*i+2]
        fx = -((cos(theta)*Adh[ind, 0]-sin(theta)*Adh[ind, 1] +
                cell[3*i] - Adh0[ind, 0]))
        fy = -((sin(theta)*Adh[ind, 0]+cos(theta)*Adh[ind, 1] +
                cell[3*i+1] - Adh0[ind, 1]))

        ad_energy += 0.5*ks(np.add(fx**2 + fy**2))

    return ad_energy

def overlap_energy(cells, Ovlaps, k_o):

    oval_energy = 0                
