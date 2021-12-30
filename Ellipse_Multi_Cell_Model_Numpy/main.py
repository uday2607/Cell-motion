import Cell_funcs
import Adh_funcs
import energy
from scipy.optimize import minimize
#from minimize_CER import minimize
import numpy as np
import plot_cells
import math

a = 6
b = 3
L = 40
Num = 1
Nadh = 64
k_plus = 0.5
dt = 1
lamda = 0.7
tau = 30
T_S = 10
k_s = 0.5
k_out_out = 50
k_in_out = k_out_out*50
k_in_in = k_in_out*50
fThreshold = 0.0005
k_b = 0.0025
k_f = 0.001
alpha = 25
a_min = a
for i in range(tau):
    a_min -= a_min*lamda*dt/tau
rng = np.random.default_rng()

# Spawn Cells and Adhesions
cells, cparams, Ovlaps = Cell_funcs.random_cells(L, a, b, Num)
#cells, cparams, Ovlaps = Cell_funcs.preset(a, b, Num)
cells0 = cells.copy()
cparams0 = cparams.copy()
Adh0, Adh = Adh_funcs.random_adhesions(a, b, cells, cparams, Nadh,
                                       k_plus, dt)

for t in range(6*(tau+tau//3+1)):
    print("Time = ", t)

    #plot in the beginning
    plot_cells.plot_ellipses(cells, cparams, Adh, Adh0, a, b, t)

    #update cells
    for num in range(Num):
        #cell extension
        if cparams[4*num+3] > 0:
            #Contract
            print('c')
            cparams[4*num:(4*num+4)], Adh[num] = Adh_funcs.contraction(cells, \
                                            num, cparams, Ovlaps, Adh, Adh0, \
                                            lamda, tau, dt, T_S, \
                                            k_out_out, k_in_out, k_in_in, k_s)
        elif cparams[4*num+3] < 0:
            #Protrusion
            print('p')
            cells[3*num:(3*num+3)], cparams[4*num:(4*num+4)], \
            Adh[num] = Adh_funcs.protrusion(cells, num, Adh, Adh0, \
                                    cparams, Ovlaps, T_S, lamda, k_s, \
                                    k_out_out, k_in_out, k_in_in, dt, tau)

        else:
            #New bonds
            print('n')
            Adh0[num], Adh[num] = Adh_funcs.one_cell_random_adh(a, b, cells,
                                    num, cparams, Adh, Adh0, Nadh,
                                    k_plus, dt)
            cparams[4*num+3] = 1

    #minimize energy
    args = (cparams, Ovlaps, Adh, Adh0,
            k_s, k_out_out, k_in_out, k_in_in)
    cells_ = cells.copy()
    cells = minimize(energy.total_energy, cells, args=args,
                    method='L-BFGS-B').x
    #cells, conv, _ = minimize(energy.total_energy, cells, args=args)
    #print(cells, cells_)

    #rotate and shift the adhesions
    cells, cparams, Adh = Adh_funcs.rotation_and_shift(cells, cells_,
                            cparams, Adh)

    # Update bonds
    for num in range(Num):
        if cparams[4*num+3] == 1:
            #first phase
            Adh[num], Adh0[num] = Adh_funcs.mature(cells[3*num:(3*num+3)],
                                    Adh[num], Adh0[num],
                                    k_s, fThreshold, dt, rng)
        elif cparams[4*num+3] != 1 and cparams[4*num+3] != 0:
            #detach
            Adh[num], Adh0[num] = Adh_funcs.detach(cells[3*num:(3*num+3)],
                                    cells0[3*num:(3*num+3)],
                                    Adh[num], Adh0[num],
                                    cparams0[4*num:(4*num+3)],
                                    k_b, k_f, alpha, a, dt, rng)
