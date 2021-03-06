from cell import *
from energy import *
from scipy.optimize import fsolve
import scipy.stats as stats
import sys

# Function to find distance
def shortest_distance(points, a, b, c):

    x1 = points[:, 0]
    y1 = points[:, 1]
    d = np.abs((a * x1 + b * y1 + c)) / (np.sqrt(a * a + b * b))
    return d

def random_adhesions(L, a, b, amin, cells, cparams,
                    Nadh, k_plus, dt, rng):

    Adh = np.zeros((cells.shape[0]//3, Nadh, 2)) - 1e8
    Adh0 = np.zeros((cells.shape[0]//3, Nadh, 2)) - 1e8

    for ind in range(cells.shape[0]//3):

        #Random numbers 
        phi = rng.normal(0, pi, Nadh)
        rho = rng.uniform(0, 1, Nadh)

        #In the system frame of reference
        x = cparams[4*ind]*np.sqrt(rho)*np.cos(phi)
        y = cparams[4*ind+1]*np.sqrt(rho)*np.sin(phi)

        #Transform into ellipse
        xp = (x*np.cos(pi/180*cparams[4*ind+2]) -
                        y*np.sin(pi/180*cparams[4*ind+2]))
        yp = (x*np.sin(pi/180*cparams[4*ind+2]) +
                        y*np.cos(pi/180*cparams[4*ind+2]))

        rands = rng.random(Nadh)

        #indices of adhesions to be formed
        adh_ind = np.arange(Nadh)[rands < k_plus*dt]

        #Store adhesions
        #From ground reference
        Adh0[ind, adh_ind, 0] = xp[adh_ind] + cells[3*ind]
        Adh0[ind, adh_ind, 1] = yp[adh_ind] + cells[3*ind+1]
        #From centre of cell(GFR)
        Adh[ind, adh_ind, 0] = xp[adh_ind]
        Adh[ind, adh_ind, 1] = yp[adh_ind]

    return Adh0, Adh

def one_cell_random_adh(a, b, cells, ind, cparams, Adh_, Adh0_,
                    Nadh, k_plus, dt, rng, flag):

    Adh = np.array(Adh_)
    Adh0 = np.array(Adh0_)

    #Indices of adhesion sites not bonded yet
    adh_ind = np.arange(Nadh)[np.all(Adh[ind]==-1e8,axis=1)]

    #Random numbers 
    phi = rng.normal(0, pi, Nadh)
    rho = rng.uniform(0, 1, Nadh)

    #In the system frame of reference
    x = cparams[4*ind]*np.sqrt(rho)*np.cos(phi)
    y = cparams[4*ind+1]*np.sqrt(rho)*np.sin(phi)

    #Transform into ellipse
    xp = (x*np.cos(pi/180*cparams[4*ind+2]) -
            y*np.sin(pi/180*cparams[4*ind+2]))
    yp = (x*np.sin(pi/180*cparams[4*ind+2]) +
            y*np.cos(pi/180*cparams[4*ind+2]))

    rands = rng.random(adh_ind.shape[0])

    #indices of adhesions to be formed
    adh_ind = adh_ind[rands < k_plus*dt]

    #Store adhesions
    #From ground reference
    Adh0[ind, adh_ind, 0] = xp[adh_ind] + cells[3*ind]
    Adh0[ind, adh_ind, 1] = yp[adh_ind] + cells[3*ind+1]
    #From centre of cell(GFR)
    Adh[ind, adh_ind, 0] = xp[adh_ind]
    Adh[ind, adh_ind, 1] = yp[adh_ind]
        
    return Adh[ind], Adh0[ind]

def rotation_and_shift(cells_, cell_p, cparams_, Adh_, Adh0_):

    cells = cells_.copy()
    cparams = cparams_.copy()
    Adh = Adh_.copy()
    Adh0 = Adh0_.copy()

    for i in range(cells.shape[0]//3):
        ind = np.arange(Adh.shape[1])[np.all(Adh[i]!=-1e8,axis=1)]

        #rotate theta
        cparams[4*i+2] += cells_[3*i+2]
        cells[3*i+2] = 0

        #rotate all the adhesions
        x, y = Adh_[i, ind, 0], Adh_[i, ind, 1]
        dtheta = pi/180*cells_[3*i+2]
        Adh[i, ind, 0] = np.cos(dtheta)*x - np.sin(dtheta)*y
        Adh[i, ind, 1] = np.sin(dtheta)*x + np.cos(dtheta)*y

    return cells, cparams, Adh, Adh0

def contraction(cells_, i, cparams_, Ovlaps_, Adh_, Adh0_, lamda, tau, dt, T_S,
                k_out_out, k_in_out, k_in_in, k_s):

    cells = np.array(cells_)
    Adh = np.array(Adh_)
    Adh0 = np.array(Adh0_)
    Ovlaps = np.array(Ovlaps_)
    cparams = np.array(cparams_)

    #Sanity check
    #dtheta should be zero
    if cells[3*i+2] != 0:
        print("Dtheta is not zero(Contraction)")
        print("Code exiting")
        exit(1)

    #compute tension due to other cells
    T = compute_tension(cells, i, cparams, Ovlaps, Adh, Adh0,
                        k_out_out, k_in_out, k_in_in, k_s)
    #Tension is less than critical tension
    if (T < T_S and T > 0):
        c = lamda*dt*(1-T/T_S)/(tau)

        # Update a and phase
        cparams[4*i] -= c*cparams[4*i]
        cparams[4*i+3] += 1

        ind = np.arange(Adh.shape[1])[np.all(Adh[i]!=-1e8,axis=1)]
        x, y = Adh[i, ind, 0], Adh[i, ind, 1]
        theta = pi/180*(cparams[4*i+2])

        Adh[i, ind, 0] = (x - c*np.cos(theta)*np.cos(theta)*x -
                          c*np.sin(theta)*np.cos(theta)*y)
        Adh[i, ind, 1] = (y - c*np.sin(theta)*np.cos(theta)*x -
                          c*np.sin(theta)*np.sin(theta)*y)
    elif (T <= 0):
        # negative tension means no stalling
        c = lamda*dt/(tau)

        # Update a and phase
        cparams[4*i] -= c*cparams[4*i]
        cparams[4*i+3] += 1

        ind = np.arange(Adh.shape[1])[np.all(Adh[i]!=-1e8,axis=1)]
        x, y = Adh[i, ind, 0], Adh[i, ind, 1]
        theta = pi/180*(cparams[4*i+2])

        Adh[i, ind, 0] = (x - c*np.cos(theta)*np.cos(theta)*x -
                         c*np.sin(theta)*np.cos(theta)*y)
        Adh[i, ind, 1] = (y - c*np.sin(theta)*np.cos(theta)*x -
                         c*np.sin(theta)*np.sin(theta)*y)
    else:
        #Tension is more than critical tension
        #No shift of adhesions

        #Update Phase
        cparams[4*i+3] += 1

    #Change phase if the contraction phase is over
    if cparams[4*i+3] > 0 and (cparams[4*i+3]-1) % tau == 0:
        cparams[4*i+3] = -1 # -ve times -> Protrusion

    return cparams[4*i:(4*i+4)], Adh[i]

def solve_center(vals, a, b, theta, h, k, xc, yc):

    x, y = vals

    if theta > np.pi:
        theta = -2*np.pi + theta

    return [(((h+xc)-x)*np.cos(theta)+((k+yc)-y)*np.sin(theta))**2/a**2+
            (((h+xc)-x)*np.sin(theta)-((k+yc)-y)*np.cos(theta))**2/b**2 - 1,
            np.arctan2(y-yc, x-xc)-theta]

def protrusion(cells, num, Adh, Adh0, cparams, Ovlaps,
               T_S, lamda, k_s, k_out_out, k_in_out,
               k_in_in, dt, tau):

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

    #velocity factor
    vf = 1.0     

    #Tension is less than critical tension
    if (T < T_S and T > 0):
        c = lamda*dt*(1-T/T_S)/(tau//10)

        # Update a and phase
        cparams[4*num] += c*cparams[4*num]
        cparams[4*num+3] -= 1

        # indices of adhesions
        ind = np.arange(Adh.shape[1])[np.logical_and(
                Adh[num, :, 1] != -1e8, Adh[num, :, 0] != -1e8)]

        #If there are FAs
        if ind.shape[0] != 0:

            # Find the perpendicular line to semi major axis
            a1 = -1/tan(cparams[4*num+2])
            b1 = -1
            c1 = (-a1*(cparams[4*num]*np.cos(cparams[4*num+2])) +
             cparams[4*num]*np.sin(cparams[4*num+2]))

            #Find the rear most adhesion
            ad = ind[np.argmax(shortest_distance(Adh[num, ind],
                    a1, b1, c1))]

            #Find the new center of the cell
            theta = cparams[4*num+2]*pi/180
            args = (cparams[4*num], cparams[4*num+1], theta,
                    Adh[num, ad, 0], Adh[num, ad, 1], cells[3*num],
                    cells[3*num+1])
            #(Solving for new center)
            x_, y_ = fsolve(solve_center, cells[3*num:3*num+2]+np.ones(2), args=args)
            dis = np.sqrt((cells[3*num] - x_)**2. + (cells[3*num+1] - y_)**2.)
            xn, yn = cells[3*num:3*num+2] + vf*10*dis/tau*np.array([np.cos(theta),
                                                                np.sin(theta)])

            x, y = (Adh[num, ind, 0], Adh[num, ind, 1])

            #Shift the adhesions inside and translate the center
            Adh[num, ind, 0] = (x + c*np.cos(theta)*np.cos(theta)*x +
                            c*np.sin(theta)*np.cos(theta)*y)
            Adh[num, ind, 1] = (y + c*np.sin(theta)*np.cos(theta)*x +
                            c*np.sin(theta)*np.sin(theta)*y)

            #Shift the center
            cells[3*num], cells[3*num+1] = xn, yn

    elif (T <= 0):
        # negative tension means no stalling
        c = lamda*dt/(tau//10)

        # Update a and phase
        cparams[4*num] += c*cparams[4*num]
        cparams[4*num+3] -= 1

        # indices of adhesions
        ind = np.arange(Adh.shape[1])[np.logical_and(
                Adh[num, :, 1] != -1e8, Adh[num, :, 0] != -1e8)]

        #If there are FAs
        if ind.shape[0] != 0:

            # Find the perpendicular line to semi major axis
            a1 = -1/tan(cparams[4*num+2])
            b1 = -1
            c1 = (-a1*(cparams[4*num]*np.cos(cparams[4*num+2])) +
             cparams[4*num]*np.sin(cparams[4*num+2]))

            #Find the rear most adhesion
            ad = ind[np.argmax(shortest_distance(Adh[num, ind],
                    a1, b1, c1))]

            #Find the new center of the cell
            theta = cparams[4*num+2]*pi/180
            args = (cparams[4*num], cparams[4*num+1], theta,
                    Adh[num, ad, 0], Adh[num, ad, 1], cells[3*num],
                    cells[3*num+1])
            #(Solving for new center)
            x_, y_ = fsolve(solve_center, cells[3*num:3*num+2]+np.ones(2), args=args)
            dis = np.sqrt((cells[3*num] - x_)**2. + (cells[3*num+1] - y_)**2.)
            xn, yn = cells[3*num:3*num+2] + vf*10*dis/tau*np.array([np.cos(theta),
                                                                np.sin(theta)])

            x, y = (Adh[num, ind, 0], Adh[num, ind, 1])

            #Shift the adhesions inside and translate the center
            Adh[num, ind, 0] = (x + c*np.cos(theta)*np.cos(theta)*x +
                            c*np.sin(theta)*np.cos(theta)*y)
            Adh[num, ind, 1] = (y + c*np.sin(theta)*np.cos(theta)*x +
                            c*np.sin(theta)*np.sin(theta)*y)

            #Shift the center
            cells[3*num], cells[3*num+1] = xn, yn
    else:
        #Tension is more than critical tension
        #No shift of adhesions

        #Update Phase
        cparams[4*num+3] -= 1

        #If there are no adhesions, the cell won't move
        #So the cell just contracts

    #Change phase if the protrusion phase is over
    if cparams[4*num+3] < 0 and abs(cparams[4*num+3]) % (tau//10+1) == 0:
        cparams[4*num+3] = 0 # +ve times -> Contraction
                             # zero time -> New Adhesions

    return cells[3*num:3*num+3], cparams[4*num:(4*num+4)], Adh[num], Adh0[num]

def mature(cell, Adh, Adh0, k_s, fTh, dt, rng):

    ind = np.arange(Adh.shape[0])[np.all(Adh !=-1e8,axis=1)]

    for i in ind:
        dis = np.sqrt(np.sum((Adh[i]+cell[:2]-Adh0[i])**2.))
        off_rate = exp(-k_s*dis/fTh)

        #detach adhesion site
        if (rng.random() < off_rate*dt):
            Adh[i] = np.array([-1e8, -1e8])
            Adh0[i] = np.array([-1e8, -1e8])

    return Adh, Adh0

def detach(cell, cell0, Adh, Adh0, cp0,
            k_b, k_f, alpha, a, dt, rng):

    ind = np.arange(Adh.shape[0])[np.all(Adh!=-1e8,axis=1)]

    for i in ind:
        x0 = ((Adh0[i, 0] - cell0[0])*np.cos(cp0[2]) +
              (Adh0[i, 1] - cell0[1])*np.sin(cp0[2]))

        #detachment rate
        k_x = k_b - (k_b - k_f)*(x0+a)/(2*a)

        dis = np.sqrt(np.sum((Adh[i]+cell[:2]-Adh0[i])**2.))
        off_rate = k_x*exp(alpha*dis/(2*a))

        #detach adhesion site
        if (rng.random() < off_rate*dt):
            Adh[i] = np.array([-1e8, -1e8])
            Adh0[i] = np.array([-1e8, -1e8])

    return Adh, Adh0
