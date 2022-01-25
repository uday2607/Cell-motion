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
    Y_MAX = 40
    Y_MIN = 0
    Num = 25
    Nadh = 100
    k_plus = 0.01
    lamda = 0.5
    tau = 30
    dt = 60/tau
    T_S = 8000
    k_s = 1.0
    k_m = 1
    k_out_out = 40.0
    k_in_out = k_out_out*1000
    k_in_in = k_in_out*1000
    fThreshold = 0.005*k_s
    k_b = 0.01
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
    #cells, cparams, Ovlaps = Cell_funcs.linear_preset(a, b, Num)
    cells, cparams, Ovlaps = Cell_funcs.custom_preset(a, b, Num)
    #cells, cparams, Ovlaps = Cell_funcs.two_linear_preset(a, b, Num)
    cparams0 = cparams.copy()
    Adh0, Adh, cAdh0, cp0 = Adh_funcs.random_adhesions(a, b, cells, cparams, Nadh,
                                           k_plus, dt)

    #CELLS = np.zeros((TIME, cells.shape[0]))
    #CPARAMS = np.zeros((TIME, cparams.shape[0]))
    #ADH = np.zeros((TIME, Adh.shape[0], Adh.shape[1], Adh.shape[2]))
    #ADH0 = np.zeros(ADH.shape)
    #CADH0 = np.zeros(ADH.shape)
    #CP0 = np.zeros(ADH.shape)

    with open("distance_data.txt", "w") as f:
        f.write("{}\t{}\t{}\t{}\n".format(cells[0], cells[1], 0.0, "n"))

    for t in range(TIME):
        print("Time = ", t)

        ## Store data
        #CELLS[t] = cells.copy()
        #CPARAMS[t] = cparams.copy()
        #ADH[t] = Adh.copy()
        #ADH0[t] = Adh0.copy()
        #CADH0[t] = cAdh0.copy()
        #CP0[t] = cp0.copy()

        #plot in the beginning
        plot_cells.plot_ellipses(cells, cparams, Adh, Adh0, a, b, Y_MAX, Y_MIN, t)

        #Find overlap between cells
        Ovlaps = Cell_funcs.find_overlaps(cells, cparams, Ovlaps)

        #update cells
        for num in range(Num):
            #cell extension
            if cparams[4*num+3] > 0:
                #Contract
                print('c')
                cparams[4*num:(4*num+4)], Adh[num], c_flag = Adh_funcs.contraction(cells, \
                                                num, cparams, Ovlaps, Adh, Adh0, \
                                                lamda, tau, dt, T_S, \
                                                k_out_out, k_in_out, k_in_in, k_s, a_min)

                #Old bonds detach
                print("detach")
                if cparams[4*num+3] == 2:
                    #first phase
                    Adh[num], Adh0[num], cAdh0[num], cp0[num] = Adh_funcs.mature(cells[3*num:(3*num+3)],
                                        Adh[num], Adh0[num], cAdh0[num], cp0[num],
                                        k_m, fThreshold, dt, rng)
                elif cparams[4*num+3] > 2:
                    #detach
                    Adh[num], Adh0[num], cAdh0[num], cp0[num] = Adh_funcs.detach(cells[3*num:(3*num+3)],
                                        cparams[4*num:4*num+4], cAdh0[num], cp0[num],
                                        Adh[num], Adh0[num],
                                        k_b, k_f, alpha, a, dt, rng)

                print("new")
                Adh0[num], Adh[num], cAdh0[num], cp0[num] = Adh_funcs.one_cell_random_adh(a, b, cells,
                                        num, cparams, Adh, Adh0, cAdh0, cp0, Nadh,
                                       0.5*k_plus, dt)

                with open("distance_data.txt", "a+") as f:
                    f.write("{}\t{}\t{}\t{}\n".format(cells[0], cells[1], t, "c"))

            elif cparams[4*num+3] < 0:
                #Protrusion
                print('p')
                cells[3*num:(3*num+3)], cparams[4*num:(4*num+4)], \
                Adh[num], p_flag = Adh_funcs.protrusion(cells, num, Adh, Adh0, cparams, Ovlaps,
                            lamda, tau, a, a_min, k_out_out, k_in_out, k_in_in, k_s, T_S)

                #Old bonds detach -> New bonds form
                print("detach")
                Adh[num], Adh0[num], cAdh0[num], cp0[num]= Adh_funcs.detach(cells[3*num:(3*num+3)],
                                    cparams[4*num:4*num+4], cAdh0[num], cp0[num],
                                    Adh[num], Adh0[num],
                                    k_b, k_f, alpha, a, dt, rng)

                #New adhesions at a smaller rate
                print("new")
                if p_flag:
                    Adh0[num], Adh[num], cAdh0[num], cp0[num] = Adh_funcs.one_cell_random_adh(a, b, cells,
                                        num, cparams, Adh, Adh0, cAdh0, cp0, Nadh,
                                        2*k_plus, dt)

                with open("distance_data.txt", "a+") as f:
                    f.write("{}\t{}\t{}\t{}\n".format(cells[0], cells[1], t, "p"))

        #minimize energy
        args = (cparams, Ovlaps, Adh, Adh0,
                k_s, k_out_out, k_in_out, k_in_in, Y_MAX, Y_MIN)
        cells_ = cells.copy()

        print(energy.total_energy(cells, cparams, Ovlaps, Adh, Adh0,
                k_s, k_out_out, k_in_out, k_in_in, Y_MAX, Y_MIN))
        soln = minimize_parallel(fun=energy.total_energy, x0=cells, args=args, jac=energy.total_energy_gradient,
        options={"maxfun" : 10**5, "maxiter" : 10**5})
        print(soln["success"])
        cells = soln.x
        print(energy.total_energy(cells, cparams, Ovlaps, Adh, Adh0,
                k_s, k_out_out, k_in_out, k_in_in, Y_MAX, Y_MIN))

        #rotate and shift the adhesions
        cells, cparams, Adh = Adh_funcs.rotation_and_shift(cells, cells_,
                                cparams, Adh)

    # Save the data
    #with open("data.npy", "wb") as f:
    #    np.save(f, CELLS)
    #    np.save(f, CPARAMS)
    #    np.save(f, ADH)
    #    np.save(f, ADH0)
    #    np.save(f, CADH0)
    #    np.save(f, CP0)
