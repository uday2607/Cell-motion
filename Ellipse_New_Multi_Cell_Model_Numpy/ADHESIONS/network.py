from Cell_funcs import *

def mature(cell, Adh, Adh0, cAdh0, cp0, k_m, fTh, dt, rng):

    ind = np.arange(Adh.shape[0])[np.all(Adh !=-1e8,axis=1)]

    for i in ind:
        if cp0[i, 2] == 0:
            dis = sqrt(np.sum((Adh[i]+cell[:2]-Adh0[i])**2.))
            off_rate = k_m*exp(-dis/fTh)

            #detach adhesion site
            p = rng.random()
            if (p < 1-off_rate*dt):
                Adh[i] = np.array([-1e8, -1e8])
                Adh0[i] = np.array([-1e8, -1e8])
                cAdh0[i] = np.array([-1e8, -1e8])
                cp0[i] = np.array([-1e8, -1e8, -1e8])
            else:
                cp0[i, 2] = 1
        elif cp0[i, 2] == 1:
            cp0[i, 2] = 2        

    return Adh, Adh0, cAdh0, cp0

def detach(cell, cparam, cAdh0, cp0, Adh, Adh0,
            k_b, k_f, alpha, a, dt, rng):

    ind = np.arange(Adh.shape[0])[np.all(Adh!=-1e8,axis=1)]

    for i in ind:
        #If adhesion is not inside, remove it
        a = cparam[0]
        b = cparam[1]
        theta = cparam[2]

        if ((((Adh[i, 0])*np.cos(theta)+(Adh[i, 1])*np.sin(theta))**2/a**2+
            ((Adh[i, 0])*np.sin(theta)-(Adh[i, 1])*np.cos(theta))**2/b**2) > 1):
            Adh[i] = np.array([-1e8, -1e8])
            Adh0[i] = np.array([-1e8, -1e8])
            cAdh0[i] = np.array([-1e8, -1e8])
            cp0[i] = np.array([-1e8, -1e8, -1e8])

        x0 = ((cAdh0[i, 0])*cos(cp0[i, 1]) +
              (cAdh0[i, 1])*sin(cp0[i, 1]))

        #detachment rate
        k_x = k_b - (k_b - k_f)*(x0+cp0[i, 0])/(2*cp0[i, 0])

        dis = sqrt(np.sum((np.around(Adh[i]+cell[:2]-Adh0[i], 4))**2.))
        #The offrate increases as cell length decreases
        off_rate = k_x*exp(alpha*dis/(2*cp0[i, 0]))

        #detach adhesion site
        p = rng.random()
        if (p < off_rate*dt) and cp0[i, 2] == 2:
            Adh[i] = np.array([-1e8, -1e8])
            Adh0[i] = np.array([-1e8, -1e8])
            cAdh0[i] = np.array([-1e8, -1e8])
            cp0[i] = np.array([-1e8, -1e8, -1e8])

    return Adh, Adh0, cAdh0, cp0
