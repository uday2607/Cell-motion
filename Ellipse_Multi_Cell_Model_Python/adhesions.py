from cell import *
from energy import *

def random_adhesions(L, a, b, cells, cparams,
                    Nadh, k_plus, dt, rng):

    Adh = np.zeros((cells.shape[0]//3, Nadh, 2)) - 1
    Adh0 = np.zeros((cells.shape[0]//3, Nadh, 2)) - 1

    for ind in range(cells.shape[0]//3):
        n = 2
        #Front and back points must be connected
        #front
        x, y = a, 0.0
        xp = (x*cos(pi/180*cparams[4*ind+2]) +
                y*sin(pi/180*cparams[4*ind+2]))
        yp = (x*sin(pi/180*cparams[4*ind+2]) +
                y*cos(pi/180*cparams[4*ind+2]))
        Adh0[ind, 0] = np.array([xp, yp]) + cells[3*ind:(3*ind+2)]
        Adh[ind, 0] = np.array([xp, yp])

        #back
        x, y = -a, 0.0
        xp = (x*cos(pi/180*cparams[4*ind+2]) +
                y*sin(pi/180*cparams[4*ind+2]))
        yp = (x*sin(pi/180*cparams[4*ind+2]) +
                y*cos(pi/180*cparams[4*ind+2]))
        Adh0[ind, 1] = np.array([xp, yp]) + cells[3*ind:(3*ind+2)]
        Adh[ind, 1] = np.array([xp, yp])

        for i in range(Nadh-2):
            flag = 0
            if (rng.random() < k_plus*dt):
                while flag == 0:
                    x, y = rng.uniform(0, L, size=2)

                    if (x**2/a**2 + y**2/b**2 <= 1):
                        xp = (x*cos(pi/180*cparams[4*ind+2]) +
                                y*sin(pi/180*cparams[4*ind+2]))
                        yp = (x*sin(pi/180*cparams[4*ind+2]) +
                                y*cos(pi/180*cparams[4*ind+2]))

                        #Store adhesions
                        Adh0[ind, n] = np.array([xp, yp]) + cells[3*ind:(3*ind+2)]
                        Adh[ind, n] = np.array([xp, yp])
                        n += 1
                        flag = 1

    return Adh0, Adh

def new_adhesion(cell, cp, Nadh, a, b, L, k_plus, dt, rng):

    Adh = np.zeros((Nadh, 2)) - 1
    Adh0 = np.zeros((Nadh, 2)) - 1

    #Front and back points must be connected
    #front
    x, y = a, 0.0
    xp = (x*cos(pi/180*cp[2]) + y*sin(pi/180*cp[2]))
    yp = (x*sin(pi/180*cp[2]) + y*cos(pi/180*cp[2]))
    Adh0[0] = np.array([xp, yp]) + cell[:2]
    Adh[0] = np.array([xp, yp])

    #back
    x, y = -a, 0.0
    xp = (x*cos(pi/180*cp[2]) + y*sin(pi/180*cp[2]))
    yp = (x*sin(pi/180*cp[2]) + y*cos(pi/180*cp[2]))
    Adh0[1] = np.array([xp, yp]) + cell[:2]
    Adh[1] = np.array([xp, yp])

    n = 2
    for i in range(Nadh-2):
        flag = 0
        if (rng.random() < k_plus*dt):
            while flag == 0:
                x, y = rng.uniform(0, L, size=2)

                if (x**2/a**2 + y**2/b**2 <= 1):
                    xp = (x*cos(pi/180*cp[2]) + y*sin(pi/180*cp[2]))
                    yp = (x*sin(pi/180*cp[2]) + y*cos(pi/180*cp[2]))

                    #Store adhesions
                    Adh0[n] = np.array([xp, yp]) + cell[:2]
                    Adh[n] = np.array([xp, yp])
                    n += 1
                    flag = 1

    return Adh, Adh0

def rotation(cells, cparams, Adh):

    for i in range(cells.shape[0]//3):
        ind = np.arange(Adh.shape[1])[np.all(Adh[i]==-1,axis=1)]
        x, y = Adh[i, ind, 0], Adh[i, ind, 1]
        dtheta = pi/180*cells[3*i+2]

        #rotate all the adhesions
        Adh[i, ind, 0] = cos(dtheta)*x - sin(dtheta)*y
        Adh[i, ind, 1] = sin(dtheta)*x + cos(dtheta)*y

        #rotate theta
        cparams[4*i+2] += cells[3*i+2]
        cells[3*i+2] = 0

    return cells, cparams, Adh

def contraction(cells, cparams, Ovlaps, Adh, lamda, tau, dt, T_S, a, b,
                k_out_out, k_in_out, k_in_in):

    for i in range(cells.shape[0]//3):
        # If no cells overlap
        if (np.all(Ovlaps[i] == -1)):
             c = lamda*dt/tau

             # Update a
             cparams[4*i] = a*exp(-c*cparams[4*i+3])
             ind = np.arange(Adh.shape[1])[np.all(Adh[i]==-1,axis=1)]
             x, y = Adh[i, ind, 0], Adh[i, ind, 1]
             theta = pi/180*cparams[4*i+2]

             Adh[i, ind, 0] = (x - c*cos(theta)*cos(theta)*x -
                                c*sin(theta)*cos(theta)*y)
             Adh[i, ind, 1] = (y - c*sin(theta)*cos(theta)*x -
                                c*sin(theta)*sin(theta)*y)
        else:
            #cells overlap

            #compute tension due to other cells
            T = compute_tension(cells, i, cparams, Ovlaps,
                                k_out_out, k_in_out, k_in_in)
            #Tension is less than critical tension
            if (T < T_S):
                c = lamda*dt*(1-T/T_S)/tau

                # Update a
                cparams[4*i] = a*exp(-c*cparams[4*i+3])

                ind = np.arange(Adh.shape[1])[np.all(Adh[i]==-1,axis=1)]
                x, y = Adh[i, ind, 0], Adh[i, ind, 1]
                theta = pi/180*cparams[4*i+2]

                Adh[i, ind, 0] = (x - c*cos(theta)*cos(theta)*x -
                                   c*sin(theta)*cos(theta)*y)
                Adh[i, ind, 1] = (y - c*sin(theta)*cos(theta)*x -
                                   c*sin(theta)*sin(theta)*y)
            elif (T <= 0):
                # negative tension means no stalling
                c = lamda*dt/tau

                # Update a
                cparams[4*i] = a*exp(-c*cparams[4*i+3])
                ind = np.arange(Adh.shape[1])[np.all(Adh[i]==-1,axis=1)]
                x, y = Adh[i, ind, 0], Adh[i, ind, 1]
                theta = pi/180*cparams[4*i+2]

                Adh[i, ind, 0] = (x - c*cos(theta)*cos(theta)*x -
                                   c*sin(theta)*cos(theta)*y)
                Adh[i, ind, 1] = (y - c*sin(theta)*cos(theta)*x -
                                   c*sin(theta)*sin(theta)*y)
            else:
                #Tension is more than critical tension
                #No shift of adhesions
                continue

    return cparams, Adh

def mature(cell, Adh, Adh0, cp, Nadh, k_s, fTh, dt, rng):

    ind = np.arange(Adh.shape[0])[np.all(Adh ==-1,axis=1)]

    for i in ind:
        dis = sqrt(np.sum((Adh[i]-Adh0[i])**2.))
        off_rate = exp(-k_s*dis/fTh)

        #detach adhesion site
        if (rng.random() < off_rate*dt):
            Adh[i] = -np.array([-1, -1])
            Adh0[i] = -np.array([-1, -1])

    return Adh, Adh0

def detach(cell, cell0, Adh, Adh0, cp, cp0,
            Nadh, k_b, k_f, k_s, alpha, dt, rng):

    ind = np.arange(Adh.shape[0])[np.all(Adh ==-1,axis=1)]

    for i in ind:
        x0 = ((Adh0[i, 0] - cell0[i])*cos(cp0[2]) +
              (Adh0[i, 1] - cell0[i+1])*sin(cp0[2]))

        #detachment rate
        k_x = k_b - (k_b - k_f)*(x0+a)/(2*a)

        dis = sqrt(np.sum((Adh[i]-Adh0[i])**2.))
        off_rate = k_x*exp(alpha*dis/a*ks/0.1)

        #detach adhesion site
        if (rng.random() < off_rate*dt):
            Adh[i] = -np.array([-1, -1])
            Adh0[i] = -np.array([-1, -1])

    return Adh, Adh0
