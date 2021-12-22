import cell
import adhesions
import energy
from scipy.optimize import minimize
import numpy as np
import plot_cells

a = 4
b = 2
L = 40
Num = 2
Nadh = 32
k_plus = 0.4
dt = 1
lamda = 0.4
tau = 30
T_S = 10
k_s = 0.1
k_out_out = 0.02
k_in_out = 1000
k_in_in = 10000
fThreshold = 0.0005
k_b = 0.02
k_f = 0.005
alpha = 25
R = 2*a
a_min = a*np.exp(-lamda*29/30)
rng = np.random.default_rng()

# Spawn Cells and Adhesions
cells, cparams, Ovlaps = cell.random_cells(L, a, b, Num, rng)
cells0 = cells.copy()
cparams0 = cparams.copy()
Adh0, Adh = adhesions.random_adhesions(L, a, b, cells, cparams, Nadh,
                                       k_plus, dt, rng)

for t in range(4*tau):

    #Find the overlap indices
    Ovlaps = cell.find_overlaps(cells, cparams, Ovlaps)

    # Cells contract
    cparams, Adh = adhesions.contraction(cells, cparams, Ovlaps, Adh,
                          lamda, tau, dt, T_S, a, b,
                          k_out_out, k_in_out, k_in_in)

    #minimize energy
    args = (cparams, Ovlaps, Adh, Adh0,
            k_s, k_out_out, k_in_out, k_in_in)
    cells = minimize(energy.total_energy, cells, args=args,
                    method='L-BFGS-B').x

    #rotate the adhesions
    cells, cparams, Adh = adhesions.rotation(cells, cparams, Adh)

    #update phase
    cparams[4*np.arange(Num)+3] += 1

    #update the bonds
    for num in range(Num):
        #cell extension
        if cparams[4*num] <= a_min:
            #reset phase
            cparams[4*num] = a
            cparams[4*num+3] = 0

            #Protrusion
            cells0[3*num:(3*num+3)] = cells[3*num:(3*num+3)].copy()
            cparams0[4*num:(4*num+4)] = cparams[4*num:(4*num+4)].copy()

            #new adhesions
            Adh[num], Adh0[num] = adhesions.new_adhesion(cells[3*num:(3*num+3)],
                                    cparams[4*num:(4*num+3)], Nadh, a, b, L, k_plus,
                                    dt, rng)
        elif cparams[4*num+3] == 1:
            #first phase
            Adh[num], Adh0[num] = adhesions.mature(cells[3*num:(3*num+3)],
                                    Adh[num], Adh0[num],
                                    cparams[4*num:(4*num+3)], Nadh,
                                    k_s, fThreshold, dt, rng)
        else:
            #detach
            Adh[num], Adh0[num] = adhesions.detach(cells[3*num:(3*num+3)],
                                    cells0[3*num:(3*num+3)],
                                    Adh[num], Adh0[num],
                                    cparams[4*num:(4*num+3)],
                                    cparams0[4*num:(4*num+3)],
                                    Nadh, k_b, k_f, k_s, alpha, a, dt, rng)
    plot_cells.plot_ellipses(cells, cparams, Adh, Adh0, a, b, t)
    print(t)
