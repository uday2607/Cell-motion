import Cell_funcs
import Adh_funcs
import energy
from scipy.optimize import minimize
from optimparallel import minimize_parallel
import numpy as np
import plot_cells
import math

if __name__ == '__main__':

    a = 6
    b = 3
    L = 40
    Num = 5
    Nadh = 64
    k_plus = 0.5
    dt = 1
    lamda = 0.5
    tau = 30
    T_S = 100
    k_s = 0.5
    k_m = 1
    k_out_out = 2
    k_in_out = k_out_out*100
    k_in_in = k_in_out*100
    fThreshold = 0.5
    k_b = 0.005
    k_f = 0.0005
    alpha = 25
    a_min = a
    for i in range(tau):
        a_min -= a_min*lamda*dt/tau    
    rng = np.random.default_rng()

    # Spawn Cells and Adhesions
    #cells, cparams, Ovlaps = Cell_funcs.random_cells(L, a, b, Num)
    #cells, cparams, Ovlaps = Cell_funcs.preset(a, b, Num)
    cells, cparams, Ovlaps = Cell_funcs.linear_preset(a, b, Num)
    cparams0 = cparams.copy()
    Adh0, Adh, cAdh0, cp0 = Adh_funcs.random_adhesions(a, b, cells, cparams, Nadh,
                                           k_plus, dt)

    for t in range(6*(tau+tau//5+1)+1):
        print("Time = ", t)

        #plot in the beginning
        plot_cells.plot_ellipses(cells, cparams, Adh, Adh0, a, b, t)

        #Find overlap between cells
        Ovlaps = Cell_funcs.find_overlaps(cells, cparams, Ovlaps)

        #update cells
        for num in range(Num):
            #cell extension
            if cparams[4*num+3] > 0:
                #Contract
                print('c')
                cparams[4*num:(4*num+4)], Adh[num] = Adh_funcs.contraction(cells, \
                                                num, cparams, Ovlaps, Adh, Adh0, \
                                                lamda, tau, dt, T_S, \
                                                k_out_out, k_in_out, k_in_in, k_s, a_min)                                          

                #Old bonds detach
                if cparams[4*num+3] == 2:
                    #first phase
                    Adh[num], Adh0[num], cAdh0[num], cp0[num] = Adh_funcs.mature(cells[3*num:(3*num+3)],
                                        Adh[num], Adh0[num], cAdh0[num], cp0[num],
                                        k_m, fThreshold/k_s, dt, rng)                  
                elif cparams[4*num+3] != 1 and cparams[4*num+3] != 0 and cparams[4*num+3] != 2:
                    #detach
                    Adh[num], Adh0[num], cAdh0[num], cp0[num] = Adh_funcs.detach(cells[3*num:(3*num+3)],
                                        cAdh0[num], cp0[num],
                                        Adh[num], Adh0[num],
                                        k_b, k_f, alpha, a, dt, rng)


            elif cparams[4*num+3] < 0:
                #Protrusion
                print('p')
                cells[3*num:(3*num+3)], cparams[4*num:(4*num+4)], \
                Adh[num] = Adh_funcs.protrusion(cells, num, Adh, Adh0, \
                                        cparams, Ovlaps, T_S, lamda, k_s, \
                                        k_out_out, k_in_out, k_in_in, dt, tau, a)

                #Old bonds detach -> New bonds form
                if cparams[4*num+3] != 0:
                    Adh[num], Adh0[num], cAdh0[num], cp0[num] = Adh_funcs.detach(cells[3*num:(3*num+3)],
                                        cAdh0[num], cp0[num],
                                        Adh[num], Adh0[num],
                                        k_b, k_f, 2*alpha, a, dt, rng)                        

                #New adhesions at a smaller rate
                Adh0[num], Adh[num], cAdh0[num], cp0[num] = Adh_funcs.one_cell_random_adh(a, b, cells,
                                        num, cparams, Adh, Adh0, cAdh0, cp0, Nadh,
                                        k_plus/10, dt, 1)                        

            else:
                #New bonds
                print('n')
                Adh0[num], Adh[num], cAdh0[num], cp0[num] = Adh_funcs.one_cell_random_adh(a, b, cells,
                                        num, cparams, Adh, Adh0, cAdh0, cp0, Nadh,
                                        k_plus, dt, 0)
                cparams[4*num+3] = 1                            

        #minimize energy
        args = (cparams, Ovlaps, Adh, Adh0,
                k_s, k_out_out, k_in_out, k_in_in)     
        cells_ = cells.copy()
        bounds = []
        eps = 1e-8
        for i in range(Num):
            ub = np.max((cells[3*i]-cparams[4*i]*np.cos(cparams[4*i+2]-eps),
                           cells[3*i]+cparams[4*i]*np.cos(cparams[4*i+2])+eps))
            lb = np.min((cells[3*i]-cparams[4*i]*np.cos(cparams[4*i+2]-eps),
                           cells[3*i]+cparams[4*i]*np.cos(cparams[4*i+2])+eps))
            bounds.append((lb, ub))
            ub = np.max((cells[3*i+1]-cparams[4*i]*np.sin(cparams[4*i+2]-eps),
                           cells[3*i+1]+cparams[4*i]*np.sin(cparams[4*i+2])+eps))
            lb = np.min((cells[3*i+1]-cparams[4*i]*np.sin(cparams[4*i+2]-eps),
                           cells[3*i+1]+cparams[4*i]*np.sin(cparams[4*i+2])+eps))               
            bounds.append((lb, ub))
            bounds.append((0, 2*np.pi))        
        bounds = tuple(bounds)

        #soln = minimize(energy.total_energy, cells, args=args,
        #                method='BFGS')
        soln = minimize_parallel(fun=energy.total_energy, x0=cells, args=args)
        cells = soln.x        

        #rotate and shift the adhesions
        cells, cparams, Adh = Adh_funcs.rotation_and_shift(cells, cells_,
                                cparams, Adh)

        # Update bonds
        #for num in range(Num):

        
