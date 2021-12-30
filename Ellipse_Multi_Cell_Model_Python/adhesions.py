from cell import *
from energy import *
from scipy.optimize import fsolve

# Function to find distance
def shortest_distance(points, a, b, c):

    x1 = points[:, 0]
    y1 = points[:, 1]
    d = np.abs((a * x1 + b * y1 + c)) / (np.sqrt(a * a + b * b))
    return d

def random_adhesions(L, a, b, cells, cparams,
                    Nadh, k_plus, dt, rng):

    Adh = np.zeros((cells.shape[0]//3, Nadh, 2)) - 1e8
    Adh0 = np.zeros((cells.shape[0]//3, Nadh, 2)) - 1e8

    for ind in range(cells.shape[0]//3):
        n = 2
        #Front and back points must be connected
        #front
        x = a
        xp = x*cos(pi/180*cparams[4*ind+2])
        yp = x*sin(pi/180*cparams[4*ind+2])
        Adh0[ind, 0] = np.array([xp, yp]) + cells[3*ind:(3*ind+2)]
        Adh[ind, 0] = np.array([xp, yp])

        #back
        x = -a
        xp = x*cos(pi/180*cparams[4*ind+2])
        yp = x*sin(pi/180*cparams[4*ind+2])
        Adh0[ind, 1] = np.array([xp, yp]) + cells[3*ind:(3*ind+2)]
        Adh[ind, 1] = np.array([xp, yp])

        for i in range(Nadh-2):

            phi = rng.normal(0, pi/2)
            rho = rng.uniform(0, 1)

            #In the system frame of reference
            x = a*sqrt(rho)*cos(phi)
            y = b*sqrt(rho)*sin(phi)

            if (rng.random() < k_plus*dt):
                #Transform into ellipse
                xp = (x*cos(pi/180*cparams[4*ind+2]) -
                                y*sin(pi/180*cparams[4*ind+2]))
                yp = (x*sin(pi/180*cparams[4*ind+2]) +
                                y*cos(pi/180*cparams[4*ind+2]))

                #Store adhesions
                #From ground reference
                Adh0[ind, i] = np.array([xp, yp]) + cells[3*ind:(3*ind+2)]
                #From centre of cell(GFR)
                Adh[ind, i] = np.array([xp, yp])

    return Adh0, Adh

def one_cell_random_adh(a, b, cells, ind, cparams, Adh_, Adh0_,
                    Nadh, k_plus, dt, rng):

    Adh = np.array(Adh_)
    Adh0 = np.array(Adh0_)

    #Indices of adhesion sites not bonded yet
    adh_ind = np.arange(Adh.shape[1])[np.all(Adh[ind]==-1e8,axis=1)]

    if adh_ind.shape[0] != 0:
        n = 2
        #Front and back points must be connected
        #front
        x = cparams[4*ind]
        xp = x*cos(cparams[4*ind+2]*pi/180)
        yp = x*sin(cparams[4*ind+2]*pi/180)
        Adh0[ind, adh_ind[0]] = np.array([xp, yp]) + cells[3*ind:(3*ind+2)]
        Adh[ind, adh_ind[0]] = np.array([xp, yp])

        #back
        x = -cparams[4*ind]
        xp = x*cos(cparams[4*ind+2]*pi/180)
        yp = x*sin(cparams[4*ind+2]*pi/180)
        Adh0[ind, adh_ind[1]] = np.array([xp, yp]) + cells[3*ind:(3*ind+2)]
        Adh[ind, adh_ind[1]] = np.array([xp, yp])

        if adh_ind.shape[0] > 2:
            for i in adh_ind[2:]:

                phi = rng.normal(0, pi/2)
                rho = rng.uniform(0, 1)

                #In the system frame of reference
                x = cparams[4*ind]*sqrt(rho)*cos(phi)
                y = cparams[4*ind+1]*sqrt(rho)*sin(phi)

                if (rng.uniform(0, 1) < k_plus*dt):
                    #Transform into ellipse
                    xp = (x*cos(cparams[4*ind+2]*pi/180) -
                                    y*sin(cparams[4*ind+2]))
                    yp = (x*sin(cparams[4*ind+2]*pi/180) +
                                    y*cos(cparams[4*ind+2]))

                    #Store adhesions
                    #From ground reference
                    Adh0[ind, i] = np.array([xp, yp]) + cells[3*ind:(3*ind+2)]
                    #From centre of cell(GFR)
                    Adh[ind, i] = np.array([xp, yp])

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
        ell = create_ellipse(cells[3*i:3*i+2], cparams[4*i:4*i+2], cparams[4*i+2])

        for n in ind:
            #Shift the adhesions
            Adh[i, n] += -cells_[3*i:3*i+2] + cell_p[3*i:3*i+2]

            #rotate all the adhesions
            x, y = Adh_[i, n, 0], Adh_[i, n, 1]
            dtheta = pi/180*cells_[3*i+2]
            Adh[i, n, 0] = cos(dtheta)*x - sin(dtheta)*y
            Adh[i, n, 1] = sin(dtheta)*x + cos(dtheta)*y

            if not ell.contains(Point(Adh[i, n, 0]+cells[3*i], Adh[i, n, 1]+cells[3*i+1])):
                Adh[i, n] = -1e8*np.ones(2)
                Adh0[i, n] = -1e8*np.ones(2)

        cells[3*i+2] = 0

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

        Adh[i, ind, 0] = (x - c*cos(theta)*cos(theta)*x -
                          c*sin(theta)*cos(theta)*y)
        Adh[i, ind, 1] = (y - c*sin(theta)*cos(theta)*x -
                          c*sin(theta)*sin(theta)*y)
    elif (T <= 0):
        # negative tension means no stalling
        c = lamda*dt/(tau)

        # Update a and phase
        cparams[4*i] -= c*cparams[4*i]
        cparams[4*i+3] += 1

        ind = np.arange(Adh.shape[1])[np.all(Adh[i]!=-1e8,axis=1)]
        x, y = Adh[i, ind, 0], Adh[i, ind, 1]
        theta = pi/180*(cparams[4*i+2])

        Adh[i, ind, 0] = (x - c*cos(theta)*cos(theta)*x -
                         c*sin(theta)*cos(theta)*y)
        Adh[i, ind, 1] = (y - c*sin(theta)*cos(theta)*x -
                         c*sin(theta)*sin(theta)*y)
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

    return [((x-h)*cos(theta)+(y-k)*sin(theta))**2/a**2+
            ((x-h)*sin(theta)-(y-k)*cos(theta))**2/b**2 - 1,
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
            ad = ind[np.argmax(shortest_distance(Adh[num, ind],
                    a1, b1, c1))]

            #Find the new center of the cell
            theta = cparams[4*num+2]*pi/180
            args = (cparams[4*num], cparams[4*num+1], theta,
                    Adh[num, ad, 0], Adh[num, ad, 1], cells[3*num],
                    cells[3*num+1])
            #(Solving for new center)
            x_, y_ = fsolve(solve_center, cells[3*num:3*num+2], args=args)
            dis = sqrt((cells[3*num] - x_)**2. + (cells[3*num+1] - y_)**2.)
            xn, yn = cells[3*num:3*num+2] + 3*dis/tau*np.array([cos(theta),
                                                                sin(theta)])

            #Shift the adhesion sites and center
            xc, yc = cells[3*num], cells[3*num+1]
            ell = create_ellipse((xc, yc), cparams[4*num:4*num+2],
                                theta*180/pi)
            for i in ind:
                x, y = (Adh[num, i, 0] + xc - xn,
                        Adh[num, i, 1] + yc - yn)

                Adh[num, i, 0] = (x + c*cos(theta)*cos(theta)*x +
                                  c*sin(theta)*cos(theta)*y)
                Adh[num, i, 1] = (y + c*sin(theta)*cos(theta)*x +
                                  c*sin(theta)*sin(theta)*y)

                #check if the adhesions are inside the shifted ellipse
                if not ell.contains(Point(Adh[num, i, 0]+xn, Adh[num, i, 1]+yn)):
                    Adh[num, i] = -1e8*np.array([1, 1])
                    Adh0[num, i] = -1e8*np.array([1, 1])

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
            ad = ind[np.argmax(shortest_distance(Adh[num, ind],
                    a1, b1, c1))]

            #Find the new center of the cell
            theta = cparams[4*num+2]*pi/180
            args = (cparams[4*num], cparams[4*num+1], theta,
                    Adh[num, ad, 0], Adh[num, ad, 1], cells[3*num],
                    cells[3*num+1])
            #(Solving for new center)
            x_, y_ = fsolve(solve_center, cells[3*num:3*num+2], args=args)
            dis = sqrt((cells[3*num] - x_)**2. + (cells[3*num+1] - y_)**2.)
            xn, yn = cells[3*num:3*num+2] + 3*dis/tau*np.array([cos(theta),
                                                                sin(theta)])

            #Shift the adhesion sites and center
            xc, yc = cells[3*num], cells[3*num+1]
            ell = create_ellipse((xc, yc), cparams[4*num:4*num+2], theta*180/pi)
            for i in ind:
                x, y = (Adh[num, i, 0] + xc - xn,
                        Adh[num, i, 1] + yc - yn)

                Adh[num, i, 0] = (x + c*cos(theta)*cos(theta)*x +
                                  c*sin(theta)*cos(theta)*y)
                Adh[num, i, 1] = (y + c*sin(theta)*cos(theta)*x +
                                  c*sin(theta)*sin(theta)*y)

                #check if the adhesions are inside the shifted ellipse
                if not ell.contains(Point(Adh[num, i, 0]+xn, Adh[num, i, 1]+yn)):
                    Adh[num, i] = -1e8*np.array([1, 1])
                    Adh0[num, i] = -1e8*np.array([1, 1])

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

        dis = sqrt(np.sum((Adh[i]+cell[:2]-Adh0[i])**2.))
        off_rate = k_x*exp(alpha*dis/(2*a))

        #detach adhesion site
        if (rng.random() < off_rate*dt):
            Adh[i] = np.array([-1e8, -1e8])
            Adh0[i] = np.array([-1e8, -1e8])

    return Adh, Adh0
