import cell
import adhesions
import energy
from scipy.optimize import minimize
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
T_S = 20
k_s = 0.5
k_out_out = 15
k_in_out = k_out_out*10
k_in_in = k_in_out*10
fThreshold = 0.0005
k_b = 0.02
k_f = 0.001
alpha = 25
a_min = a
for i in range(tau):
    a_min -= a_min*lamda*dt/tau
rng = np.random.default_rng()

# Spawn Cells and Adhesions
cells, cparams, Ovlaps = cell.random_cells(L, a, b, Num, rng)
#cells, cparams, Ovlaps = cell.preset(L, a, b, Num)
cells0 = cells.copy()
cparams0 = cparams.copy()
Adh0, Adh = adhesions.random_adhesions(L, a, b, a_min, cells, cparams, Nadh,
                                       k_plus, dt, rng)

for t in range(6*(tau+tau//10+1)+1):
    print("Time = ", t)

    #plot in the beginning
    plot_cells.plot_ellipses(cells, cparams, Adh, Adh0, a, b, t)

    #update cells
    for num in range(Num):
        #cell extension
        if cparams[4*num+3] > 0:
            #Contract
            print('c')
            cparams[4*num:(4*num+4)], Adh[num] = adhesions.contraction(cells, \
                                            num, cparams, Ovlaps, Adh, Adh0, \
                                            lamda, tau, dt, T_S, \
                                            k_out_out, k_in_out, k_in_in, k_s)

            #New adhesions at much smaller rate
            '''Adh[num], Adh0[num] = adhesions.one_cell_random_adh(a, b, cells,
                                    num, cparams, Adh, Adh0, Nadh,
                                    k_plus/(5*cparams[4*num+3]+1), dt, rng)'''
        elif cparams[4*num+3] < 0:
            #Protrusion
            print('p')
            cells[3*num:(3*num+3)], cparams[4*num:(4*num+4)], \
            Adh[num], Adh0[num] = adhesions.protrusion(cells, num, Adh, Adh0, \
                                    cparams, Ovlaps, T_S, lamda, k_s, \
                                    k_out_out, k_in_out, k_in_in, dt, tau)

            #New adhesions at much smaller rate
            Adh[num], Adh0[num] = adhesions.one_cell_random_adh(a, b, cells,
                                    num, cparams, Adh, Adh0, Nadh,
                                    k_plus/5, dt, rng, 1)
        else:
            #New bonds
            print('n')
            Adh[num], Adh0[num] = adhesions.one_cell_random_adh(a, b, cells,
                                    num, cparams, Adh, Adh0, Nadh,
                                    k_plus, dt, rng, 0)
            cparams[4*num+3] = 1
            cells0 = cells.copy()      

        if cparams[4*num+3] == 2:
            #first phase
            Adh[num], Adh0[num] = adhesions.mature(cells[3*num:(3*num+3)],
                                    Adh[num], Adh0[num],
                                    k_s, fThreshold, dt, rng)
        elif cparams[4*num+3] != 1 and cparams[4*num+3] != 2 and cparams[4*num+3] != 0:
            #detach
            Adh[num], Adh0[num] = adhesions.detach(cells[3*num:(3*num+3)],
                                    cells0[3*num:(3*num+3)],
                                    Adh[num], Adh0[num],
                                    cparams0[4*num:(4*num+3)],
                                    k_b, k_f, alpha, a, dt, rng)    

    #minimize energy
    args = (cparams, Ovlaps, Adh, Adh0,
            k_s, k_out_out, k_in_out, k_in_in)
    cells_ = cells.copy()
    cells = minimize(energy.total_energy, cells+np.ones(3*Num), args=args,
                    method='L-BFGS-B', jac='3-point').x        

    #rotate and shift the adhesions
    cells, cparams, Adh, Adh0 = adhesions.rotation_and_shift(cells, cells_,
                            cparams, Adh, Adh0)                         
    if Num == 1:
        print(np.arctan2(cells[1] - cells_[1], cells[0] - cells_[0])*180/np.pi
        , cparams[3])                                   

    # Update bonds
    #for num in range(Num):
        
