from Cell_funcs import *
from energy import compute_tension

"""Contraction of cell -> contraction of Adh sites"""
@nb.jit(nopython = True, nogil = True)
def contraction(cells_, num, cparams_, Ovlaps, Adh_, Adh0_, lamda, tau,
                dt, T_S, E_const, k_s, a_min):

    #copy the arrays
    cells = cells_.copy()
    cparams = cparams_.copy()
    Adh = Adh_.copy()
    Adh0_ = Adh0_.copy()

    #Sanity check
    #dtheta should be zero
    if np.any(cells[3*np.arange(cells.shape[0]//3)+2] != 0):
        print("Dtheta is not zero(Contraction)")
        while 1:
            print("Stop the code")

    #compute tension due to other cells pulling outwards
    T = compute_tension(cells_, num, cparams_, Ovlaps, Adh_, Adh0_,
                        E_const, k_s)

    #Tension is less than critical tension
    if (T < T_S and T > 0):
        print("Tension:", T)
        c = lamda*dt*(1-T/T_S)/(tau)

        # Update a and phase
        cparams[4*num] -= c*cparams[4*num]
        cparams[4*num+3] += 1
        c_flag = True

        # indices of adhesions
        ind = np.arange(Adh.shape[1])[np.logical_and(
                Adh[num, :, 1] != -1e8, Adh[num, :, 0] != -1e8)]
        theta = cparams[4*num+2]

        for i in ind:
            x, y = Adh[num, i, 0], Adh[num, i, 1]

            Adh[num, i, 0] = (x - c*cos(theta)*cos(theta)*x -
                              c*sin(theta)*cos(theta)*y)                 
            Adh[num, i, 1] = (y - c*sin(theta)*cos(theta)*x -
                              c*sin(theta)*sin(theta)*y)         

    elif (T <= 0):
        print("Tension:", T)
        # negative tension means no stalling
        c = lamda*dt/(tau)

        # Update a and phase
        cparams[4*num] -= c*cparams[4*num]
        cparams[4*num+3] += 1
        c_flag = True

        # indices of adhesions
        ind = np.arange(Adh.shape[1])[np.logical_and(
                Adh[num, :, 1] != -1e8, Adh[num, :, 0] != -1e8)]
        theta = cparams[4*num+2]

        for i in ind:
            x, y = Adh[num, i, 0], Adh[num, i, 1]

            Adh[num, i, 0] = (x - c*cos(theta)*cos(theta)*x -
                              c*sin(theta)*cos(theta)*y)
            Adh[num, i, 1] = (y - c*sin(theta)*cos(theta)*x -
                              c*sin(theta)*sin(theta)*y)
    else:
        #Tension is more than critical tension
        #No shift of adhesions

        #Update Phase
        cparams[4*num+3] += 1
        c_flag = False
        print("Tension (critical):", T)

    #Change phase if the contraction phase is over
    if cparams[4*num+3] > 0 and cparams[4*num] < a_min:
        cparams[4*num+3] = -1 # -ve times -> Protrusion

    return cparams[4*num:(4*num+4)], Adh[num], c_flag
""""""