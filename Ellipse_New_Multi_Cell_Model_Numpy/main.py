import Cell_funcs
import Adh_funcs
import energy
from scipy.optimize import minimize
from optimparallel import minimize_parallel
import numpy as np
import plot_cells
import os
import math

if __name__ == '__main__':

    a = 4
    b = 2
    L = 1000
    Num = 10
    Nadh = 64
    k_plus = 0.025
    lamda = 0.2
    tau = 30
    dt = 60/tau
    T_S_c = 20000
    T_S_p = 20000
    k_s = 50.0
    k_m = 1
    E_const = np.zeros((3, 3))
    E_const[2, 2] = -0.025
    E_const[1, 1] = 2.0
    E_const[1, 0] = 5000.0
    E_const[0, 1] = 5000.0
    E_const[0, 0] = 10000.0
    fThreshold = 40.0/k_s
    k_b = 0.025
    k_f = 0.0025
    alpha = 25
    a_min = a
    for i in range(tau):
        a_min -= a_min*lamda*dt/tau
    rng = np.random.default_rng()
    TIME = 600

    # To save the data of the simulations

    # Spawn Cells and Adhesions
    #cells, cparams, Ovlaps = Cell_funcs.random_cells(L, a, b, Num)
    cells, cparams, Ovlaps = Cell_funcs.linear_preset(a, b, Num)
    #cells, cparams, Ovlaps = Cell_funcs.two_linear_preset(a, b, Num)
    #cells, cparams, Ovlaps = Cell_funcs.four_linear_preset(a, b, Num)
    cparams0 = cparams.copy()
    Adh0, Adh, cAdh0, cp0 = Adh_funcs.random_adhesions(cells, cparams, Nadh,
                                           10*k_plus, dt, rng)

    CELLS = np.zeros((TIME, cells.shape[0]))
    CPARAMS = np.zeros((TIME, cparams.shape[0]))
    ADH = np.zeros((TIME, Adh.shape[0], Adh.shape[1], Adh.shape[2]))
    ADH0 = np.zeros(ADH.shape)
    CADH0 = np.zeros(ADH.shape)
    CP0 = np.zeros((ADH.shape[0], ADH.shape[1], ADH.shape[2], ADH.shape[3]+1))

    with open("distance_data.txt", "w") as f:
        f.write("{}\t{}\t{}\t{}\n".format(cells[0], cells[1], 0.0, "n"))

    for t in range(TIME):
        print("Time = ", t)

        #plot in the beginning
        plot_cells.plot_ellipses(cells, cparams, Adh, Adh0, a, b, t) 

        ## Store data
        CELLS[t] = cells.copy()
        CPARAMS[t] = cparams.copy()
        ADH[t] = Adh.copy()
        ADH0[t] = Adh0.copy()
        CADH0[t] = cAdh0.copy()
        CP0[t] = cp0.copy()

        #Find overlap between cells
        Ovlaps = Cell_funcs.find_overlaps(cells, cparams, Ovlaps)

        #Make array copies
        cparams_ = cparams.copy()
        Adh_ = Adh.copy()
        Adh0_ = Adh0.copy()
        cAdh0_ = cAdh0.copy()
        cp0_ = cp0.copy()

        #update cells
        for num in range(Num):
            #cell extension
            if cparams_[4*num+3] > 0:
                #Contract
                print("c")
                cparams[4*num:(4*num+4)], Adh[num], c_flag = Adh_funcs.contraction(cells, \
                                                num, cparams_, Ovlaps, Adh_, Adh0_, \
                                                lamda, tau, dt, T_S_c, \
                                                E_const, k_s, a_min)                              

                #Old bonds detach
                if cparams_[4*num+3] > 1:
                    #first phase
                    Adh[num], Adh0[num], cAdh0[num], cp0[num] = Adh_funcs.mature(cells[3*num:(3*num+3)],
                                        Adh[num], Adh0[num], cAdh0[num], cp0[num],
                                        k_m, fThreshold, dt, rng)
                    #detach
                    Adh[num], Adh0[num], cAdh0[num], cp0[num] = Adh_funcs.detach(cells[3*num:(3*num+3)],
                                        cparams[4*num:4*num+4], cAdh0[num], cp0[num],
                                        Adh[num], Adh0[num],
                                        k_b, k_f, alpha, a, dt, rng)                    

                #New adhesions at a smaller rate
                Adh0[num], Adh[num], cAdh0[num], cp0[num] = Adh_funcs.one_cell_random_adh(cells,
                                        num, cparams, Adh, Adh0, cAdh0, cp0, Nadh,
                                       k_plus, dt, 0, rng)

                with open("distance_data.txt", "a+") as f:
                    f.write("{}\t{}\t{}\t{}\n".format(cells[0], cells[1], t, "c"))

            elif cparams_[4*num+3] < 0:
                #Protrusion
                print('p')
                cells[3*num:(3*num+3)], cparams[4*num:(4*num+4)], \
                p_flag = Adh_funcs.protrusion(cells, num, Adh_, Adh0_, cparams_, Ovlaps,
                            dt, tau, a, a_min, E_const, k_s, T_S_p)

                #New bonds form
                if p_flag:
                    Adh0[num], Adh[num], cAdh0[num], cp0[num] = Adh_funcs.one_cell_random_adh(cells,
                                        num, cparams, Adh, Adh0, cAdh0, cp0, Nadh,
                                        10*k_plus, dt, 1, rng)

                with open("distance_data.txt", "a+") as f:
                    f.write("{}\t{}\t{}\t{}\n".format(cells[0], cells[1], t, "p"))

        #minimize energy
        args = (cparams, Ovlaps, Adh, Adh0,
                k_s, E_const)
        cells_ = cells.copy()

        #plot_cells.plot_ellipses(cells, cparams, Adh, Adh0, a, b, 2*t)
        print("Energy (before): ", energy.total_energy(cells, cparams, Ovlaps, Adh, Adh0,
                k_s, E_const))
               
        soln = minimize(fun=energy.total_energy, x0=cells, args=args, jac=energy.total_energy_gradient,
        options={"maxiter" : 10**5}, method="L-BFGS-B")
        #soln = minimize(fun=energy.total_energy, x0=cells, args=args,
        #options={"maxiter" : 10**5}, method="Nelder-Mead")
        print("Convergence: ", soln["success"])
        if soln["success"] == False:
            print(soln["message"])
            #os.system("say "+soln["message"])
        cells = soln.x
        print("Energy (after): ", energy.total_energy(cells, cparams, Ovlaps, Adh, Adh0,
                k_s, E_const))    

        #rotate and shift the adhesions
        cells, cparams, Adh = Adh_funcs.rotation_and_shift(cells, cells_,
                                cparams, Adh)             

    # Save the data
    with open("data.npy", "wb") as f:
        np.save(f, CELLS)
        np.save(f, CPARAMS)
        np.save(f, ADH)
        np.save(f, ADH0)
        np.save(f, CADH0)
        np.save(f, CP0)
