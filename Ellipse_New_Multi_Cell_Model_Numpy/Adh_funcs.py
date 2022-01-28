from Cell_funcs import *
from energy import *
#From: https://github.com/Nicholaswogan/NumbaMinpack
#from NumbaMinpack import lmdif, minpack_sig
#from numba import cfunc
from scipy.optimize import fsolve

@nb.jit(nopython = True, nogil = True)
def in_ellipse(x, y, xc, yc, a, b, theta):

    if (((x-xc)*cos(theta)+(y-yc)*sin(theta))**2/a**2 +
        ((x-xc)*sin(theta)-(y-yc)*cos(theta))**2/b**2) < 1:
        return 1
    else:
        return 0


"""New Adhesions"""
@nb.jit(nopython = True, nogil = True)
def random_adhesions(a, b, cells, cparams,
                    Nadh, k_plus, dt):

    Adh = np.zeros((cells.shape[0]//3, Nadh, 2)) - 1e8
    Adh0 = np.zeros((cells.shape[0]//3, Nadh, 2)) - 1e8
    cAdh0 = np.zeros((cells.shape[0]//3, Nadh, 2)) - 1e8
    cp0 = np.zeros((cells.shape[0]//3, Nadh, 3)) - 1e8

    for ind in range(cells.shape[0]//3):
        rands = uniform_double(0, 1, Nadh)

        #indices of adhesions to be formed
        adh_ind = np.arange(Nadh)[rands < k_plus*dt]

        for ad in adh_ind:

            # Find a point from the grid
            found = 0

            while found == 0:

                phi = uniform_double(0, 2*np.pi, 1)[0]
                rho = uniform_double(0.25, 1.0, 1)[0]
                x_ = cparams[4*ind]*sqrt(rho)*cos(phi)
                y_ = cparams[4*ind+1]*sqrt(rho)*sin(phi)
                x = x_*cos(cparams[4*ind+2]) - y_*sin(cparams[4*ind+2])
                y = x_*sin(cparams[4*ind+2]) + y_*cos(cparams[4*ind+2])
                x = round(4*(x_+cells[3*ind]))/4
                y = round(4*(y_+cells[3*ind+1]))/4

                if (in_ellipse(x, y, cells[3*ind], cells[3*ind+1],
                    cparams[4*ind], cparams[4*ind+1], cparams[4*ind+2]) and
                   not in_ellipse(x, y, cells[3*ind], cells[3*ind+1],
                    cparams[4*ind]/2, cparams[4*ind+1]/2, cparams[4*ind+2])):

                    # Check if there is an adhesion point from this site
                    if not np.any(np.logical_and(Adh[ind, :, 0] == x,
                            Adh[ind, :, 1] == y)):
                        found = 1

            #Store adhesions
            #From ground reference
            Adh0[ind, ad, 0] = x
            Adh0[ind, ad, 1] = y
            Adh[ind, ad, 0] = x - cells[3*ind]
            Adh[ind, ad, 1] = y - cells[3*ind+1]
            cAdh0[ind, ad] = Adh[ind, ad].copy()
            cp0[ind, ad, 0] = cparams[4*ind]
            cp0[ind, ad, 1] = cparams[4*ind+2]
            cp0[ind, ad, 2] = 0

    return Adh0, Adh, cAdh0, cp0

@nb.jit(nopython = True, nogil = True)
def one_cell_random_adh(a, b, cells, ind, cparams, Adh, Adh0, cAdh0, cp0,
                    Nadh, k_plus, dt):

    #Indices of adhesion sites not bonded yet
    adh_ind = np.arange(Adh.shape[1])[np.logical_and(
        Adh[ind, :, 0] == -1e8, Adh[ind, :, 1] == -1e8)]

    rands = uniform_double(0, 1, adh_ind.shape[0])

    #indices of adhesions to be formed
    adh_ind = adh_ind[rands < k_plus*dt]

    for ad in adh_ind:

        # Find a point from the grid
        found = 0

        while found == 0:
            phi = uniform_double(0, 2*np.pi, 1)[0]
            rho = uniform_double(0.25, 1.0, 1)[0]
            x_ = cparams[4*ind]*sqrt(rho)*cos(phi)
            y_ = cparams[4*ind+1]*sqrt(rho)*sin(phi)
            x = x_*cos(cparams[4*ind+2]) - y_*sin(cparams[4*ind+2])
            y = x_*sin(cparams[4*ind+2]) + y_*cos(cparams[4*ind+2])
            x = round(4*(x_+cells[3*ind]))/4
            y = round(4*(y_+cells[3*ind+1]))/4

            if (in_ellipse(x, y, cells[3*ind], cells[3*ind+1],
                cparams[4*ind], cparams[4*ind+1], cparams[4*ind+2]) and
               not in_ellipse(x, y, cells[3*ind], cells[3*ind+1],
                cparams[4*ind]/2, cparams[4*ind+1]/2, cparams[4*ind+2])):

                # Check if there is an adhesion point from this site
                if not np.any(np.logical_and(Adh[ind, :, 0] == x,
                        Adh[ind, :, 1] == y)):
                    found = 1

        #Store adhesions
        #From ground reference
        Adh0[ind, ad, 0] = x
        Adh0[ind, ad, 1] = y
        Adh[ind, ad, 0] = x - cells[3*ind]
        Adh[ind, ad, 1] = y - cells[3*ind+1]
        cAdh0[ind, ad] = Adh[ind, ad].copy()
        cp0[ind, ad, 0] = cparams[4*ind]
        cp0[ind, ad, 1] = cparams[4*ind+2]
        cp0[ind, ad, 2] = 0

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
@nb.jit(nopython = True, nogil = True)
def contraction(cells_, num, cparams_, Ovlaps, Adh_, Adh0_, lamda, tau,
                dt, T_S, k_out_out, k_in_out, k_in_in, k_s, a_min):

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
                        k_out_out, k_in_out, k_in_in, k_s)

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
@nb.jit(nopython = True, nogil = True)
def protrusion(cells_, num, Adh_, Adh0_, cparams_, Ovlaps, lamda, tau, a, a_min,
            k_out_out, k_in_out, k_in_in, k_s, T_S):

    #copy the arrays
    cells = cells_.copy()
    cparams = cparams_.copy()
    Adh = Adh_.copy()
    Adh0 = Adh0_.copy()

    #Sanity check
    #dtheta should be zero
    if np.any(cells[3*np.arange(cells.shape[0]//3)+2] != 0):
        print("Dtheta is not zero(Contraction)")
        while 1:
            print("Stop the code")

    #compute tension due to other cells pulling outwards
    T = -1.0*compute_tension(cells_, num, cparams_, Ovlaps, Adh_, Adh0_,
                        k_out_out, k_in_out, k_in_in, k_s)
    #(-ve sign for Tension to account for opposite direction of gradient)

    if (T < T_S and T > 0):
        print("Tension:", T)
        #Adhesion speed
        c = 10/tau*(1 - T/T_S)

        # Update a and phase
        cparams[4*num] += c*(a - a_min)
        cparams[4*num+3] -= 1
        p_flag = True

        # indices of adhesions
        ind = np.arange(Adh.shape[1])[np.logical_and(
            Adh[num, :, 1] != -1e8, Adh[num, :, 0] != -1e8)]

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
            theta = cparams[4*num+2]

            #Find the new center of the cell
            x_ad = Adh[num, ad, 0]*cos(theta) + Adh[num, ad, 1]*sin(theta)

            #(Solving for new center)
            xn = cells[3*num] + c*(a+x_ad)*cos(theta)
            yn = cells[3*num+1] + c*(a+x_ad)*sin(theta)

            #Shift the adhesion sites and center
            for i in ind:
                x, y = (Adh[num, i, 0], Adh[num, i, 1])

                #Shift the adhesions inside and translate the center
                Adh[num, i, 0] = x - xn + cells[3*num]
                Adh[num, i, 1] = y - yn + cells[3*num+1]

            cells[3*num], cells[3*num+1] = xn, yn

    elif T < 0.0:
        print("Tension:", T)
        #Adhesion speed
        c = 10/tau

        # Update a and phase
        cparams[4*num] += c*(a - a_min)
        cparams[4*num+3] -= 1
        p_flag = True

        # indices of adhesions
        ind = np.arange(Adh.shape[1])[np.logical_and(
            Adh[num, :, 1] != -1e8, Adh[num, :, 0] != -1e8)]

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
            theta = cparams[4*num+2]

            #Find the new center of the cell
            x_ad = Adh[num, ad, 0]*cos(theta) + Adh[num, ad, 1]*sin(theta)

            #(Solving for new center)
            xn = cells[3*num] + c*(a+x_ad)*cos(theta)
            yn = cells[3*num+1] + c*(a+x_ad)*sin(theta)

            #Shift the adhesion sites and center
            for i in ind:
                x, y = (Adh[num, i, 0], Adh[num, i, 1])

                #Shift the adhesions inside and translate the center
                Adh[num, i, 0] = x - xn + cells[3*num]
                Adh[num, i, 1] = y - yn + cells[3*num+1]

            cells[3*num], cells[3*num+1] = xn, yn

    else:
        print("Tension (critical) :", T)
        cparams[4*num+3] -= 1
        p_flag = False

    #Change phase if the protrusion phase is over
    if cparams[4*num+3] < 0 and cparams[4*num] >= a:
        cparams[4*num+3] = 1 # +ve times -> Contraction

    return cells[3*num:3*num+3], cparams[4*num:(4*num+4)], Adh[num], p_flag

def mature(cell, Adh, Adh0, cAdh0, cp0, k_m, fTh, dt, rng):

    ind = np.arange(Adh.shape[0])[np.all(Adh !=-1e8,axis=1)]
    eps = 1e-8

    for i in ind:
        if cp0[i, 2] == 0:
            dis = sqrt(np.sum((Adh[i]+cell[:2]-Adh0[i])**2.))
            off_rate = k_m*exp(-dis/fTh)

            #detach adhesion site
            if (rng.random() < 1-off_rate*dt):
                Adh[i] = np.array([-1e8, -1e8])
                Adh0[i] = np.array([-1e8, -1e8])
                cAdh0[i] = np.array([-1e8, -1e8])
                cp0[i] = np.array([-1e8, -1e8, -1e8])
            else:
                cp0[i, 2] = 1

    return Adh, Adh0, cAdh0, cp0

def detach(cell, cparam, cAdh0, cp0, Adh, Adh0,
            k_b, k_f, alpha, a, dt, rng):

    ind = np.arange(Adh.shape[0])[np.all(Adh!=-1e8,axis=1)]

    for i in ind:
        #If adhesion is not inside, remove it
        a = cparam[0]
        b = cparam[1]
        theta = cparam[2]

        if (((((Adh[i, 0])*np.cos(theta)+(Adh[i, 1])*np.sin(theta))**2/a**2+
            ((Adh[i, 0])*np.sin(theta)-(Adh[i, 1])*np.cos(theta))**2/b**2) > 1) or
           (((Adh[i, 0])*np.cos(theta)+(Adh[i, 1])*np.sin(theta))**2/(0.5*a)**2+
            ((Adh[i, 0])*np.sin(theta)-(Adh[i, 1])*np.cos(theta))**2/(0.5*b)**2) < 1):

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
        if (rng.random() < off_rate*dt):
            Adh[i] = np.array([-1e8, -1e8])
            Adh0[i] = np.array([-1e8, -1e8])
            cAdh0[i] = np.array([-1e8, -1e8])
            cp0[i] = np.array([-1e8, -1e8, -1e8])
        else:
            #update if they don't detach
            cp0[i, 2] = 1

    return Adh, Adh0, cAdh0, cp0
