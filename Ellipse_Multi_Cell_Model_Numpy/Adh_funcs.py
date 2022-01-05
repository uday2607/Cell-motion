from Cell_funcs import *
from energy import *
#From: https://github.com/Nicholaswogan/NumbaMinpack
from NumbaMinpack import lmdif, minpack_sig
from numba import cfunc

"""New Adhesions"""
@nb.jit(nopython = True, nogil = True)
def random_adhesions(a, b, cells, cparams,
                    Nadh, k_plus, dt):

    Adh = np.zeros((cells.shape[0]//3, Nadh, 2)) - 1e8
    Adh0 = np.zeros((cells.shape[0]//3, Nadh, 2)) - 1e8

    for ind in range(cells.shape[0]//3):
        n = 2
        #Front and back points must be connected
        #front
        x = a
        xp = x*cos(cparams[4*ind+2])
        yp = x*sin(cparams[4*ind+2])
        Adh0[ind, 0] = np.array([xp, yp]) + cells[3*ind:(3*ind+2)]
        Adh[ind, 0] = np.array([xp, yp])

        #back
        x = -a
        xp = x*cos(cparams[4*ind+2])
        yp = x*sin(cparams[4*ind+2])
        Adh0[ind, 1] = np.array([xp, yp]) + cells[3*ind:(3*ind+2)]
        Adh[ind, 1] = np.array([xp, yp])

        for i in range(Nadh-2):

            phi = uniform_double(0, float(2*pi), 1)[0]
            rho = uniform_double(0, 1, 1)[0]

            #In the system frame of reference
            x = a*sqrt(rho)*cos(phi)
            y = b*sqrt(rho)*sin(phi)

            if (uniform_double(0, 1, 1)[0] < k_plus*dt):
                #Transform into ellipse
                xp = (x*cos(cparams[4*ind+2]) -
                                y*sin(cparams[4*ind+2]))
                yp = (x*sin(cparams[4*ind+2]) +
                                y*cos(cparams[4*ind+2]))

                #Store adhesions
                #From ground reference
                Adh0[ind, i] = np.array([xp, yp]) + cells[3*ind:(3*ind+2)]
                #From centre of cell(GFR)
                Adh[ind, i] = np.array([xp, yp])

    return Adh0, Adh

@nb.jit(nopython = True, nogil = True)
def one_cell_random_adh(a, b, cells, ind, cparams, Adh, Adh0,
                    Nadh, k_plus, dt):

    #Indices of adhesion sites not bonded yet
    adh_ind = np.arange(Adh.shape[1])[np.logical_and(
        Adh[ind, :, 0] == -1e8, Adh[ind, :, 1] == -1e8)]

    if adh_ind.shape[0] != 0:
        n = 2
        #Front and back points must be connected
        #front
        x = a
        xp = x*cos(cparams[4*ind+2])
        yp = x*sin(cparams[4*ind+2])
        Adh0[ind, adh_ind[0]] = np.array([xp, yp]) + cells[3*ind:(3*ind+2)]
        Adh[ind, adh_ind[0]] = np.array([xp, yp])

        #back
        x = -a
        xp = x*cos(cparams[4*ind+2])
        yp = x*sin(cparams[4*ind+2])
        Adh0[ind, adh_ind[1]] = np.array([xp, yp]) + cells[3*ind:(3*ind+2)]
        Adh[ind, adh_ind[1]] = np.array([xp, yp])

        if adh_ind.shape[0] > 2:
            for i in adh_ind:

                phi = uniform_double(0, float(2*pi), 1)[0]
                rho = uniform_double(0, 1, 1)[0]

                #In the system frame of reference
                x = a*sqrt(rho)*cos(phi)
                y = b*sqrt(rho)*sin(phi)

                if (uniform_double(0, 1, 1)[0] < k_plus*dt):
                    #Transform into ellipse
                    xp = (x*cos(cparams[4*ind+2]) -
                                    y*sin(cparams[4*ind+2]))
                    yp = (x*sin(cparams[4*ind+2]) +
                                    y*cos(cparams[4*ind+2]))

                    #Store adhesions
                    #From ground reference
                    Adh0[ind, i] = np.array([xp, yp]) + cells[3*ind:(3*ind+2)]
                    #From centre of cell(GFR)
                    Adh[ind, i] = np.array([xp, yp])

    return Adh0[ind], Adh[ind]
""""""

"""Rotation and shift of Adhesion sites"""
@nb.jit(nopython = True, nogil = True)
def rotation_and_shift(cells, cell_p, cparams, Adh):

    for i in range(cells.shape[0]//3):
        ind = np.arange(Adh.shape[1])[np.logical_and(
            Adh[i, :, 0] != -1e8, Adh[i, :, 1] != -1e8)]

        if ind.shape[0] != 0:
            for n in ind:

                #rotate all the adhesions
                x, y = Adh[i, n, 0], Adh[i, n, 1]
                dtheta = cells[3*i+2]
                Adh[i, n, 0] = np.cos(dtheta)*x - np.sin(dtheta)*y
                Adh[i, n, 1] = np.sin(dtheta)*x + np.cos(dtheta)*y

                #Shift the adhesions
                #Adh[i, n, 0] += -(cells[3*i] - cell_p[3*i])
                #Adh[i, n, 1] += -(cells[3*i+1] - cell_p[3*i+1])

        #rotate theta
        cparams[4*i+2] += cells[3*i+2]
        cells[3*i+2] = 0

    return cells, cparams, Adh
""""""

"""Contraction of cell -> contraction of Adh sites"""
@nb.jit(nopython = True, nogil = True)
def contraction(cells, num, cparams, Ovlaps, Adh, Adh0, lamda, tau,
                dt, T_S, k_out_out, k_in_out, k_in_in, k_s):

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

    #Change phase if the contraction phase is over
    if cparams[4*num+3] > 0 and (cparams[4*num+3]-1) % tau == 0:
        cparams[4*num+3] = -1 # -ve times -> Protrusion

    return cparams[4*num:(4*num+4)], Adh[num]
""""""

"""Protrusion of the cell -> Complicated dynamics"""
#Helper functions
@cfunc(minpack_sig)
def solve_center(vals, fvec, args):

    #Function arguments
    #0, 1, 2, 3, 4, 5, 6
    #a, b, theta, h, k, xc, yc

    #Variables
    #0, 1
    #x, y

    #Adjust theta so that it's in arctan2 range
    if args[2] > np.pi:
        args[2] = -2*np.pi + args[2]

    #func vals
    fvec[0] = (((vals[0]-args[3])*cos(args[2])+
                (vals[1]-args[4])*sin(args[2]))**2/args[0]**2+
              ((vals[0]-args[3])*sin(args[2])-
                (vals[1]-args[4])*cos(args[2]))**2/args[1]**2 - 1)
    fvec[1] = np.arctan2(vals[1]-args[6], vals[0]-args[5])-args[2]

# address in memory to solve_center
solve_center_ptr = solve_center.address
neqs = 2

# Function to find shortest distance from a line
@nb.jit(nopython = True, nogil = True)
def shortest_distance(points, a, b, c):

    x1 = points[:, 0]
    y1 = points[:, 1]
    d = np.abs((a * x1 + b * y1 + c)) / (np.sqrt(a * a + b * b))

    return d

#Protrusion function
@nb.jit(nopython = True, nogil = True)
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

    # Find the perpendicular line to semi major axis
    a1 = -1/tan(cparams[4*num+2])
    b1 = -1
    c1 = (-a1*(cparams[4*num]*cos(cparams[4*num+2])) +
         cparams[4*num]*sin(cparams[4*num+2]))

    #Tension is less than critical tension
    if (T < T_S and T > 0):
        c = lamda*dt*(1-T/T_S)/(tau//3)

        # Update a and phase
        cparams[4*num] += c*cparams[4*num]
        cparams[4*num+3] -= 1

        # indices of adhesions
        ind = np.arange(Adh.shape[1])[np.logical_and(
                Adh[num, :, 1] != -1e8, Adh[num, :, 0] != -1e8)]

        #If there are FAs
        if ind.shape[0] != 0:

            #Find the rear most adhesion
            adh_c = Adh[num]
            ad = ind[np.argmax(shortest_distance(adh_c[ind],
                    a1, b1, c1))]

            #Find the new center of the cell
            theta = cparams[4*num+2]
            args = np.array((cparams[4*num], cparams[4*num+1], theta,
                    Adh[num, ad, 0], Adh[num, ad, 1], cells[3*num],
                    cells[3*num+1]))
            xsol, fvec, success, info = lmdif(solve_center_ptr,
                                    cells[3*num:3*num+2], neqs, args)
            #(Solving for new center)
            x_, y_ = xsol
            dis = sqrt((cells[3*num] - x_)**2. + (cells[3*num+1] - y_)**2.)
            print(dis, fvec, success, info)
            xn, yn = cells[3*num:3*num+2] + 3*dis/tau*np.array([cos(theta),
                                                                sin(theta)])

            #Shift the adhesion sites and center
            for i in ind:
                x, y = (Adh[num, i, 0], Adh[num, i, 1])

                #Shift the adhesions inside and translate the center
                Adh[num, i, 0] = (x + c*cos(theta)*cos(theta)*x +
                                  c*sin(theta)*cos(theta)*y)
                Adh[num, i, 1] = (y + c*sin(theta)*cos(theta)*x +
                                  c*sin(theta)*sin(theta)*y)

                #Adh[num, i, 0], Adh[num, i, 1] = (Adh[num, i, 0] + xc - xn,
                #        Adh[num, i, 1] + yc - yn)

            cells[3*num], cells[3*num+1] = xn, yn

    elif (T <= 0):
        # negative tension means no stalling
        c = lamda*dt/(tau//3)

        # Update a and phase
        cparams[4*num] += c*cparams[4*num]
        cparams[4*num+3] -= 1

        # indices of adhesions
        ind = np.arange(Adh.shape[1])[np.logical_and(
                Adh[num, :, 1] != -1e8, Adh[num, :, 0] != -1e8)]

        #If there are FAs
        if ind.shape[0] != 0:

            #Find the rear most adhesion
            adh_c = Adh[num]
            ad = ind[np.argmax(shortest_distance(adh_c[ind],
                    a1, b1, c1))]

            #Find the new center of the cell
            theta = cparams[4*num+2]
            args = np.array((cparams[4*num], cparams[4*num+1], theta,
                    Adh[num, ad, 0], Adh[num, ad, 1], cells[3*num],
                    cells[3*num+1]))
            xsol, fvec, success, info = lmdif(solve_center_ptr,
                                    cells[3*num:3*num+2], neqs, args)
            #(Solving for new center)
            x_, y_ = xsol
            dis = sqrt((cells[3*num] - x_)**2. + (cells[3*num+1] - y_)**2.)
            print(dis, fvec, success, info)
            xn, yn = cells[3*num:3*num+2] + 3*dis/tau*np.array([cos(theta),
                                                                sin(theta)])

            #Shift the adhesion sites and center
            for i in ind:
                x, y = (Adh[num, i, 0], Adh[num, i, 1])

                #Shift the adhesions inside and translate the center
                Adh[num, i, 0] = (x + c*cos(theta)*cos(theta)*x +
                                  c*sin(theta)*cos(theta)*y)
                Adh[num, i, 1] = (y + c*sin(theta)*cos(theta)*x +
                                  c*sin(theta)*sin(theta)*y)

                #Adh[num, i, 0], Adh[num, i, 1] = (Adh[num, i, 0] + xc - xn,
                #        Adh[num, i, 1] + yc - yn)

            cells[3*num], cells[3*num+1] = xn, yn
    else:
        #Tension is more than critical tension
        #No shift of adhesions

        #Update Phase
        cparams[4*num+3] -= 1

        #If there are no adhesions, the cell won't move

    #Change phase if the protrusion phase is over
    if cparams[4*num+3] < 0 and abs(cparams[4*num+3]) % (tau//3+1) == 0:
        cparams[4*num+3] = 0 # non -ve times -> Contraction

    return cells[3*num:3*num+3], cparams[4*num:(4*num+4)], Adh[num]

def mature(cell, Adh, Adh0, k_s, fTh, dt, rng):

    ind = np.arange(Adh.shape[0])[np.all(Adh !=-1e8,axis=1)]

    for i in ind:
        dis = sqrt(np.sum((Adh[i]+cell[:2]-Adh0[i])**2.))
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
        x0 = ((Adh0[i, 0] - cell0[0])*cos(cp0[2]) +
              (Adh0[i, 1] - cell0[1])*sin(cp0[2]))

        #detachment rate
        k_x = k_b - (k_b - k_f)*(x0+a)/(2*a)

        dis = sqrt(np.sum((np.around(Adh[i]+cell[:2]-Adh0[i], 4))**2.))
        try:
            off_rate = k_x*exp(alpha*dis/(2*a))
        except:
            print(Adh[i], cell[:2], Adh0[i])

        #detach adhesion site
        if (rng.random() < off_rate*dt):
            Adh[i] = np.array([-1e8, -1e8])
            Adh0[i] = np.array([-1e8, -1e8])

    return Adh, Adh0
