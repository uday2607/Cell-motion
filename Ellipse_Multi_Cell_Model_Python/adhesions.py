from cell import *

def random_adhesions(L, a, b, cells, Nadh, k_plus, dt, rng):

    Adh = np.zeros((cells.shape[0]//3, Nadh, 2)) - 1
    Adh0 = np.zeros((cells.shape[0]//3, Nadh, 2)) - 1

    for ind in range(cells.shape[0]//3):
        n = 0
        for i in range(Nadh):
            flag = 0
            if (rng.random() < k_plus*dt):
                while flag == 0:
                    x, y = rng.uniform(0, L, size=2)

                    if (x**2/a**2 + y**2/b**2 <= 1):
                        xp = x*cos(pi/180*cell[3*ind+2]) + y*sin(pi/180*cell[3*ind+2])
                        yp = x*sin(pi/180*cell[3*ind+2]) + y*cos(pi/180*cell[3*ind+2])

                        #Store adhesions
                        Adh0[ind, n] = np.array([xp, yp]) + cell[3*ind:3*ind+2]
                        Adh[ind, n] = np.array([xp, yp])
                        n += 1
                        flag = 1

    return Adh0, Adh

def contraction(cells, cparams, Ovlaps, Adh, lamda, tau, dt, T_S, a, b):

    for i in range(cells.shape[0]//3):
        # If no cells overlap
        if (np.all(Ovlaps[i] == -1)):
             c = lamda*cparams[3*i+2]*dt/tau

             # Update a
             cparams[3*i] = a*exp(-c)

             for ind in np.arange(Adh.shape[1])[Adh[i] != np.array([-1, -1])]:
                 x, y = Adh[i, ind]
                 theta = pi/180*cell[3*i+2]

                 Adh[i, ind, 0] = (x - c*cos(theta)*cos(theta)*x -
                                    c*sin(theta)*cos(theta)*y)
                 Adh[i, ind, 1] = (y - c*sin(theta)*cos(theta)*x -
                                    c*sin(theta)*sin(theta)*y)
        else:
            #cells overlap

            #compute tension due to other cells
            T = compute_tension(cell, cparam[i], Ovlaps[i], cells, cparams)
            #Tension is less than critical tension
            if (T < T_S):
                c = lamda*cparams[3*i+2]*dt*(1-T/T_S)/tau

                # Update a
                cparams[3*i] = a*exp(-c)

                for ind in np.arange(Adh.shape[1])[Adh[i] != np.array([-1, -1])]:
                    x, y = Adh[i, ind]
                    theta = pi/180*cell[3*i+2]

                    Adh[i, ind, 0] = (x - c*cos(theta)*cos(theta)*x -
                                       c*sin(theta)*cos(theta)*y)
                    Adh[i, ind, 1] = (y - c*sin(theta)*cos(theta)*x -
                                       c*sin(theta)*sin(theta)*y)
            else:
                #Tension is more than critical tension
                #No shift of adhesions
                continue
