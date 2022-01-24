from Cell_funcs import *
from energy import *
#From: https://github.com/Nicholaswogan/NumbaMinpack
#from NumbaMinpack import lmdif, minpack_sig
#from numba import cfunc
from scipy.optimize import fsolve

"""New Adhesions"""
@nb.jit(nopython = True, nogil = True)
def random_adhesions(a, b, cells, cparams,
                    Nadh, k_plus, dt):

    Adh = np.zeros((cells.shape[0]//3, Nadh, 2)) - 1e8
    Adh0 = np.zeros((cells.shape[0]//3, Nadh, 2)) - 1e8
    cAdh0 = np.zeros((cells.shape[0]//3, Nadh, 2)) - 1e8
    cp0 = np.zeros((cells.shape[0]//3, Nadh, 2)) - 1e8

    for ind in range(cells.shape[0]//3):

        #Random numbers 
        phi = uniform_double(0, 2*pi, Nadh)
        rho = uniform_double(0.5, 1, Nadh)

        #In the system frame of reference
        x = cparams[4*ind]*np.sqrt(rho)*np.cos(phi)
        y = cparams[4*ind+1]*np.sqrt(rho)*np.sin(phi)

        #Transform into ellipse
        xp = (x*np.cos(cparams[4*ind+2]) -
                        y*np.sin(cparams[4*ind+2]))
        yp = (x*np.sin(cparams[4*ind+2]) +
                        y*np.cos(cparams[4*ind+2]))

        rands = uniform_double(0, 1, Nadh)

        #indices of adhesions to be formed
        adh_ind = np.arange(Nadh)[rands < k_plus*dt]

        #Store adhesions
        #From ground reference
        Adh0_ = Adh0[ind].copy()
        Adh_ = Adh[ind].copy()
        cAdh0_ = cAdh0[ind].copy()
        cp0_ = cp0[ind].copy()
        Adh0_[adh_ind, 0] = xp[adh_ind] + cells[3*ind]
        Adh0_[adh_ind, 1] = yp[adh_ind] + cells[3*ind+1]
        #From centre of cell(GFR)
        Adh_[adh_ind, 0] = xp[adh_ind]
        Adh_[adh_ind, 1] = yp[adh_ind]
        cAdh0_[adh_ind] = Adh_[adh_ind]
        cp0_[adh_ind, 0] = cparams[4*ind]
        cp0_[adh_ind, 1] = cparams[4*ind+2]
        Adh[ind]= Adh_ 
        Adh0[ind] = Adh0_
        cAdh0[ind] = cAdh0_
        cp0[ind] = cp0_

    return Adh0, Adh, cAdh0, cp0

@nb.jit(nopython = True, nogil = True)
def one_cell_random_adh(a, b, cells, ind, cparams, Adh, Adh0, cAdh0, cp0,
                    Nadh, k_plus, dt, flag):

    #Indices of adhesion sites not bonded yet
    adh_ind = np.arange(Adh.shape[1])[np.logical_and(
        Adh[ind, :, 0] == -1e8, Adh[ind, :, 1] == -1e8)]

    #Random numbers 
    if flag == 2:
        phi = uniform_double(0, 2*pi, Nadh)
        rho = uniform_double(0.5, 1, Nadh) 
    elif flag == 1:
        phi = uniform_double(0, 2*pi, Nadh)
        rho = uniform_double(0.5, 1, Nadh)   
    else:
        phi = uniform_double(0, 2*pi, Nadh)
        rho = uniform_double(0.5, 1, Nadh)

    #In the system frame of reference
    x = cparams[4*ind]*np.sqrt(rho)*np.cos(phi)
    y = cparams[4*ind+1]*np.sqrt(rho)*np.sin(phi)

    #Transform into ellipse
    xp = (x*np.cos(cparams[4*ind+2]) -
            y*np.sin(cparams[4*ind+2]))
    yp = (x*np.sin(cparams[4*ind+2]) +
            y*np.cos(cparams[4*ind+2]))

    rands = uniform_double(0, 1, adh_ind.shape[0])

    #indices of adhesions to be formed
    adh_ind = adh_ind[rands < k_plus*dt]

    #Store adhesions
    #From ground reference
    Adh0_ = Adh0[ind].copy()
    Adh_ = Adh[ind].copy()
    cAdh0_ = cAdh0[ind].copy()
    cp0_ = cp0[ind].copy()
    Adh0_[adh_ind, 0] = xp[adh_ind] + cells[3*ind]
    Adh0_[adh_ind, 1] = yp[adh_ind] + cells[3*ind+1]
    #From centre of cell(GFR)
    Adh_[adh_ind, 0] = xp[adh_ind]
    Adh_[adh_ind, 1] = yp[adh_ind]
    cAdh0_[adh_ind] = Adh_[adh_ind] 
    cp0_[adh_ind, 0] = cparams[4*ind]
    cp0_[adh_ind, 1] = cparams[4*ind+2]
    Adh[ind]= Adh_ 
    Adh0[ind] = Adh0_
    cAdh0[ind] = cAdh0_
    cp0[ind] = cp0_
    
    return Adh0[ind], Adh[ind], cAdh0[ind], cp0[ind]

""""""

"""Rotation and shift of Adhesion sites"""
@nb.jit(nopython = True, nogil = True)
def rotation_and_shift(cells, cell_p, cparams, Adh):

    for i in range(cells.shape[0]//3):
        ind = np.arange(Adh.shape[1])[np.logical_and(
            Adh[i, :, 0] != -1e8, Adh[i, :, 1] != -1e8)]

        if ind.shape[0] != 0:

            #rotate all the adhesions
            Adh_ = Adh[i].copy()
            x, y = Adh_[ind, 0], Adh_[ind, 1]
            dtheta = cells[3*i+2]
            Adh_[ind, 0] = np.cos(dtheta)*x - np.sin(dtheta)*y
            Adh_[ind, 1] = np.sin(dtheta)*x + np.cos(dtheta)*y
            Adh[i] = Adh_

        #rotate theta
        cparams[4*i+2] += cells[3*i+2]
        cells[3*i+2] = 0

    return cells, cparams, Adh
""""""

"""Contraction of cell -> contraction of Adh sites"""
#@nb.jit(nopython = True, nogil = True)
def contraction(cells, num, cparams, Ovlaps, Adh, Adh0, lamda, tau,
                dt, T_S, k_out_out, k_in_out, k_in_in, k_s, a_min):

    #Sanity check
    #dtheta should be zero
    if np.any(cells[3*np.arange(cells.shape[0]//3)+2] != 0):
        print("Dtheta is not zero(Contraction)")
        while 1:
            print("Stop the code")

    #compute tension due to other cells pulling outwards
    T = compute_tension(cells, num, cparams, Ovlaps, Adh, Adh0,
                        k_out_out, k_in_out, k_in_in, k_s)

    #Tension is less than critical tension
    if (T < T_S and T > 0):
        c = lamda*dt*(1-T/T_S)/(tau)

        # Update a and phase
        cparams[4*num] -= c*cparams[4*num]
        cparams[4*num+3] += 1

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
        # negative tension means no stalling
        c = lamda*dt/(tau)

        # Update a and phase
        cparams[4*num] -= c*cparams[4*num]
        cparams[4*num+3] += 1

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
        print("Tension:", T)

    #Change phase if the contraction phase is over
    if cparams[4*num+3] > 0 and cparams[4*num] < a_min:
        cparams[4*num+3] = -1 # -ve times -> Protrusion

    return cparams[4*num:(4*num+4)], Adh[num]
""""""

"""Protrusion of the cell -> Complicated dynamics"""
#Helper functions
@nb.jit(nopython = True, nogil = True)
def solve_center(vals, a, b, theta, h, k, xc, yc):

    x, y = vals

    if theta > np.pi:
        theta = -2*np.pi + theta

    return [(((h+xc)-x)*np.cos(theta)+((k+yc)-y)*np.sin(theta))**2/a**2+
            (((h+xc)-x)*np.sin(theta)-((k+yc)-y)*np.cos(theta))**2/b**2 - 1,
            np.arctan2(y-yc, x-xc)-theta]

# Function to find shortest distance from a line
@nb.jit(nopython = True, nogil = True)
def shortest_distance(points, a, b, c):

    x1 = points[:, 0]
    y1 = points[:, 1]
    d = np.abs((a * x1 + b * y1 + c)) / (np.sqrt(a * a + b * b))

    return d

#Protrusion function
#@nb.jit(nopython = True, nogil = True)
def protrusion(cells, num, Adh, Adh0, cparams, Ovlaps,
               T_S, lamda, k_s, k_out_out, k_in_out,
               k_in_in, dt, tau, a):

    #Sanity check
    #dtheta should be zero
    if np.any(cells[3*np.arange(cells.shape[0]//3)+2] != 0):
        print("Dtheta is not zero(Contraction)")
        while 1:
            print("Stop the code")

    #compute tension due to other cells pushing inwards
    T = -1*compute_tension(cells, num, cparams, Ovlaps, Adh, Adh0,
                        k_out_out, k_in_out, k_in_in, k_s)
    #(-ve sign for Tension to account for opposite direction of gradient)

    #Adhesion speed
    vf = 1.0

    #Tension is less than critical tension
    if (T < T_S and T > 0):
        c = lamda*(1-lamda)*dt*(1-T/T_S)/(tau//5)
        p_flag = 1

        # Update a and phase
        cparams[4*num] += c*cparams[4*num]
        cparams[4*num+3] -= 1

        # indices of adhesions
        if p_flag:
            ind = np.arange(Adh.shape[1])[np.logical_and(
                Adh[num, :, 1] != -1e8, Adh[num, :, 0] != -1e8)]
        else:
            ind = np.empty((0, 2))        

        #If there are FAs
        if ind.shape[0] != 0:

            # Find the perpendicular line to semi major axis
            a1 = -1/tan(cparams[4*num+2]+1e-8)
            b1 = -1
            c1 = (-a1*(cparams[4*num]*cos(cparams[4*num+2])) +
                cparams[4*num]*sin(cparams[4*num+2]))

            #Find the rear most adhesion
            adh_c = Adh[num]
            ad = ind[np.argmax(shortest_distance(adh_c[ind],
                    a1, b1, c1))]

            #Find the new center of the cell
            theta = cparams[4*num+2]
            args = (cparams[4*num], cparams[4*num+1], theta,
                    Adh[num, ad, 0], Adh[num, ad, 1], cells[3*num],
                    cells[3*num+1])
            xsol = fsolve(solve_center, cells[3*num:3*num+2]+np.ones(2), args=args)
            #(Solving for new center)
            x_, y_ = xsol
            dis = sqrt((cells[3*num] - x_)**2. + (cells[3*num+1] - y_)**2.)
            xn, yn = cells[3*num:3*num+2] + vf*5*dis/tau*np.array([cos(theta),
                                                                sin(theta)])

            #Shift the adhesion sites and center
            for i in ind:
                x, y = (Adh[num, i, 0], Adh[num, i, 1])

                #Shift the adhesions inside and translate the center
                Adh[num, i, 0] = (x + c*cos(theta)*cos(theta)*x +
                                  c*sin(theta)*cos(theta)*y)
                Adh[num, i, 1] = (y + c*sin(theta)*cos(theta)*x +
                                  c*sin(theta)*sin(theta)*y)

            cells[3*num], cells[3*num+1] = xn, yn

    elif (T <= 0):
        # negative tension means no stalling
        c = lamda*dt/(tau//5)
        p_flag = 1

       # Update a and phase
        cparams[4*num] += c*cparams[4*num]
        cparams[4*num+3] -= 1

        # indices of adhesions
        if p_flag:
            ind = np.arange(Adh.shape[1])[np.logical_and(
                Adh[num, :, 1] != -1e8, Adh[num, :, 0] != -1e8)]
        else:
            ind = np.empty((0, 2))

        #If there are FAs
        if ind.shape[0] != 0:

            # Find the perpendicular line to semi major axis
            a1 = -1/tan(cparams[4*num+2]+1e-8)
            b1 = -1
            c1 = (-a1*(cparams[4*num]*cos(cparams[4*num+2])) +
                cparams[4*num]*sin(cparams[4*num+2]))

            #Find the rear most adhesion
            adh_c = Adh[num]
            ad = ind[np.argmax(shortest_distance(adh_c[ind],
                    a1, b1, c1))]

            #Find the new center of the cell
            theta = cparams[4*num+2]
            args = (cparams[4*num], cparams[4*num+1], theta,
                    Adh[num, ad, 0], Adh[num, ad, 1], cells[3*num],
                    cells[3*num+1])
            xsol = fsolve(solve_center, cells[3*num:3*num+2]+np.ones(2), args=args)
            #(Solving for new center)
            x_, y_ = xsol
            dis = sqrt((cells[3*num] - x_)**2. + (cells[3*num+1] - y_)**2.)
            xn, yn = cells[3*num:3*num+2] + vf*5*dis/tau*np.array([cos(theta),
                                                                sin(theta)])

            #Shift the adhesion sites and center
            for i in ind:
                x, y = (Adh[num, i, 0], Adh[num, i, 1])

                #Shift the adhesions inside and translate the center
                Adh[num, i, 0] = (x + c*cos(theta)*cos(theta)*x +
                                  c*sin(theta)*cos(theta)*y)
                Adh[num, i, 1] = (y + c*sin(theta)*cos(theta)*x +
                                  c*sin(theta)*sin(theta)*y)

            cells[3*num], cells[3*num+1] = xn, yn
    else:
        #Tension is more than critical tension
        #No shift of adhesions

        #Update Phase
        cparams[4*num+3] -= 1
        print("Tension:", T)

        #If there are no adhesions, the cell won't move

    #Change phase if the protrusion phase is over
    if cparams[4*num+3] < 0 and cparams[4*num] > a:
        cparams[4*num+3] = 0 # non -ve times -> Contraction

    return cells[3*num:3*num+3], cparams[4*num:(4*num+4)], Adh[num]

def mature(cell, Adh, Adh0, cAdh0, cp0, k_m, fTh, dt, rng):

    ind = np.arange(Adh.shape[0])[np.all(Adh !=-1e8,axis=1)]
    eps = 1e-8

    for i in ind:
        dis = sqrt(np.sum((Adh[i]+cell[:2]-Adh0[i])**2.))
        off_rate = k_m*exp(-dis/fTh)

        #detach adhesion site
        if (rng.random() < 1-off_rate*dt):
            Adh[i] = np.array([-1e8, -1e8])
            Adh0[i] = np.array([-1e8, -1e8])
            cAdh0[i] = np.array([-1e8, -1e8])
            cp0[i] = np.array([-1e8, -1e8])

    return Adh, Adh0, cAdh0, cp0

def detach(cell, cAdh0, cp0, Adh, Adh0,
            k_b, k_f, alpha, a, dt, rng):

    ind = np.arange(Adh.shape[0])[np.all(Adh!=-1e8,axis=1)]

    for i in ind:
        x0 = ((cAdh0[i, 0])*cos(cp0[i, 1]) +
              (cAdh0[i, 1])*sin(cp0[i, 1]))

        #detachment rate
        k_x = k_b - (k_b - k_f)*(x0+cp0[i, 0])/(2*cp0[i, 0])

        dis = sqrt(np.sum((np.around(Adh[i]+cell[:2]-Adh0[i], 4))**2.))
        #The offrate increases as cell length decreases
        off_rate = k_x*exp(alpha*dis/(2*cp0[i, 0]))

        #detach adhesion site
        if (rng.random() < off_rate*dt):
            Adh[i] = np.array([-1e8, -1e8])
            Adh0[i] = np.array([-1e8, -1e8])
            cAdh0[i] = np.array([-1e8, -1e8])
            cp0[i] = np.array([-1e8, -1e8])

    return Adh, Adh0, cAdh0, cp0
