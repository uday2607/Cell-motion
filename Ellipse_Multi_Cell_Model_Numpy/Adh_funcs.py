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

            phi = uniform_double(0, float(2*pi))[0]
            rho = uniform_double(0, 1)[0]

            #In the system frame of reference
            x = a*sqrt(rho)*cos(phi)
            y = b*sqrt(rho)*sin(phi)

            if (uniform_double(0, 1)[0] < k_plus*dt):
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
""""""

"""Rotation and shift of Adhesion sites"""
@nb.jit(nopython = True, nogil = True)
def rotation_and_shift(cells, cell_p, cparams, Adh):

    for i in range(cells.shape[0]//3):
        ind = np.arange(Adh.shape[1])[np.logical_and(
            Adh[0, :] != -1e8, Adh[1, :] != -1e8)]

        #Shift the adhesions 
        Adh[i, ind] += -(cells[3*i:3*i+2] - cell_p[3*i:3*i+2])

        #rotate all the adhesions
        x, y = Adh_[i, ind, 0], Adh_[i, ind, 1]
        dtheta = cells[3*i+2]
        Adh[i, ind, 0] = np.cos(dtheta)*x - np.sin(dtheta)*y
        Adh[i, ind, 1] = np.sin(dtheta)*x + np.cos(dtheta)*y

        #rotate theta
        cparams[4*i+2] += cells[3*i+2]
        cells[3*i+2] = 0

    return cells, cparams, Adh    

"""Contraction of cell -> contraction of Adh sites"""
@nb.jit(nopython = True, nogil = True)
def contraction(cells, cparams, Ovlaps, Adh, Adh0, lamda, tau, 
                dt, T_S, k_out_out, k_in_out, k_in_in, k_s):

    #Sanity check
    #dtheta should be zero
    if np.any(cells[3*np.arange(cells.shape[0]//3)+2] != 0):
        print("Dtheta is not zero(Contraction)")
        while 1:
            print("Stop the code")

    for i in range(cells.shape[0]//3):
        '''if cparams[4*i+3] == -1:
            print("Cell {} is not contracted".format(i))
            continue'''
        
        #compute tension due to other cells
        T = compute_tension(cells, i, cparams, Ovlaps, Adh, Adh0,
                            k_out_out, k_in_out, k_in_in, k_s)
        #Tension is less than critical tension
        if (T < T_S and T > 0):
            c = lamda*dt*(1-T/T_S)/(tau)

            # Update a and phase
            cparams[4*i] -= c*cparams[4*i]
            cparams[4*i+3] += 1

            ind = np.arange(Adh.shape[1])[np.logical_and(
                    Adh[0, :] != -1e8, Adh[1, :] != -1e8)]
            x, y = Adh[i, ind, 0], Adh[i, ind, 1]
            theta = cparams[4*i+2]+cells[3*i+2]

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

            ind = np.arange(Adh.shape[1])[np.logical_and(
                    Adh[0, :] != -1e8, Adh[1, :] != -1e8)]
            x, y = Adh[i, ind, 0], Adh[i, ind, 1]
            theta = cparams[4*i+2]+cells[3*i+2]

            Adh[i, ind, 0] = (x - c*cos(theta)*cos(theta)*x -
                             c*sin(theta)*cos(theta)*y)
            Adh[i, ind, 1] = (y - c*sin(theta)*cos(theta)*x -
                             c*sin(theta)*sin(theta)*y)
        else:
            #Tension is more than critical tension
            #No shift of adhesions
            continue

    return cparams, Adh
""""""    
"""Protrusion of the cell -> Complicated dynamics"""

#Helper functions
@cfunc(minpack_sig)
def solve_center(vals, fvec, args):

    #Function arguments
    a, b, theta, h, k, xc, yc = args
    #Variables
    x, y = vals 

    #Adjust theta so that it's in arctan2 range
    if theta > np.pi:
        theta = -2*np.pi + theta

    #func vals
    fvec[0] = (((x-h)*cos(theta)+(y-k)*sin(theta))**2/a**2+
            ((x-h)*sin(theta)-(y-k)*cos(theta))**2/b**2 - 1)
    fvec[1] = np.arctan2(y-yc, x-xc)-theta)

# address in memory to solve_center
solve_center_ptr = solve_center.address

# Function to find shortest distance from a line
@nb.jit(nopython = True, nogil = True)
def shortest_distance(points, a, b, c):

    x1 = points[:, 0]
    y1 = points[:, 1]
    d = np.abs((a * x1 + b * y1 + c)) / (np.sqrt(a * a + b * b))
    
    return d

#Protrusion function
@nb.jit(nopython = True, nogil = True)
def protrusion(cells, num, Adh, Adh0, cparams, 
               Nadh, a, b, k_plus, k_out_out, k_in_out, 
               k_in_in, dt, tau, rng):

    # Make copies of array
    cells_ = np.array(cells)
    cparams_ = np.array(cparams)
    cell = np.array(cells[3*num:3*num+3])
    cp = np.array(cparams[4*num:4*num+4])          

    # indices of adhesions
    ind = np.arange(Adh.shape[0])[np.all(Adh !=-1e8,axis=1)]

    # If there are FAs
    if ind.size != 0:

        # Find the perpendicular line to semi major axis
        a1 = -1/tan(cp[2])
        b1 = -1
        c1 = (-a1*(cp[0]*cos(cp[2])) +
             cp[0]*sin(cp[2]))

        #adhesion on the left most side
        ad = ind[np.argmax(shortest_distance(Adh[ind],
                    a1, b1, c1))]

        #Calculate the new position of the center
        #(The left most point should intersect the cell)
        theta = cp[2]*pi/180
        h = Adh[ad][0] + cell[0]
        k = Adh[ad][1] + cell[1]
        xc, yc = cell[:2]
        aval, bval = a, b
        args = (aval, bval, theta, h, k, xc, yc)
        xc, yc = fsolve(solve_center, cell[:2], args=args)
        cells_[3*num:3*num+2] = np.array([xc, yc])
        cparams_[4*num] = aval

        # Find overlap indices
        Ovlaps = np.zeros((Adh.shape[0], Adh.shape[0]))
        Ovlaps = find_overlaps(cells, cparams, Ovlaps)

        # Check if protruded cell overlaps      
        if (overlap_energy(cells_, cparams_, num, Ovlaps,
        k_out_out, k_in_out, k_in_in)[0] > 0):

            #Find the correct value of 'a'
            a_max = a
            a_min = cp[0]
            aval = (a_max + a_min)*0.5

            #Iterate till a value is found
            while (abs(a_max - a_min) > 1e-8):
                xc, yc = cell[:2]
                args = (aval, bval, theta, h, k, xc, yc)
                xc, yc = fsolve(solve_center, cell[:2], args=args)
                cells_[3*num:3*num+2] = np.array([xc, yc])
                cparams_[4*num] = aval

                # Find overlap indices
                Ovlaps = np.zeros((Adh.shape[0], Adh.shape[0]))
                Ovlaps = find_overlaps(cells, cparams, Ovlaps)

                #check if this cell overlaps
                if (overlap_energy(cells_, cparams_, num, Ovlaps,
                    k_out_out, k_in_out, k_in_in)[0] > 0):
                    a_max = aval 
                    a_min = a_min
                    aval = (a_max + a_min)*0.5
                    print("Right", abs(a_max - a_min))
                else:
                    a_max = a_max
                    a_min = aval
                    aval = (a_max + a_min)*0.5
                    print("Left", abs(a_max - a_min))

            #check if the cell remains in the contracted phase
            if (abs(aval - cp[0]) > 1e-8):
                #contraction phase
                cp[3] = 0
            elif cp[3] != -1:
                #no contraction phase
                cp[3] = -1
            else:
                #Don't update 'a'                        
                return cell, cp, Adh, Adh0
                
        else:
            #Update the coordinates
            cell[:2] = np.array([xc, yc])
            cp[3] = 0
            cp[0] = aval

        #Translate all adhesion sites to new cell
        h = aval - cparams[4*num]
        cp, Adh = virtual_extension(cp, Adh, h)    

        #Check if any adhesions are inside
        ell = create_ellipse(cell[:2], cp[:2], theta*180/pi)
        for i in ind:
            #inside the ellipse
            xad, yad = Adh[i] + cells[3*num:3*num+2] - cell[:2]
            if ell.contains(Point(xad, yad)):
                Adh[i] = np.array([xad, yad])
            else:
                Adh[i] = -np.array([1e8, 1e8]) 
                Adh0[i] = -np.array([1e8, 1e8])                        

        #reevaluate how many adhesions are not present in the cell
        indn = np.arange(Adh.shape[0])[np.all(Adh ==-1e8,axis=1)]

        #Check if any other adhesion is possible
        if (((indn.shape != 0) and (cparams[4*num+3] == 0))):
            #create a adhesion in the front of the cell
            x, y = aval, 0.0
            xp = (x*cos(cp[2]) - y*sin(cp[2]))
            yp = (x*sin(cp[2]) + y*cos(cp[2]))
            Adh0[indn[0]] = np.array([xp, yp]) + cell[:2]
            Adh[indn[0]] = np.array([xp, yp])

            #create a adhesion in the back of the cell
            x, y = -aval, 0.0
            xp = (x*cos(cp[2]) - y*sin(cp[2]))
            yp = (x*sin(cp[2]) + y*cos(cp[2]))
            Adh0[indn[1]] = np.array([xp, yp]) + cell[:2]
            Adh[indn[1]] = np.array([xp, yp])

            #populate other adhesions
            for ad in indn[2:]:
                phi = rng.uniform(0, float(2*pi))
                rho = rng.uniform(0, 1)

                #In the system frame of reference
                x = cp[0]*sqrt(rho)*cos(phi)
                y = b*sqrt(rho)*sin(phi)

                if (rng.random() < k_plus*dt):
                    #Transform into ellipse
                    xp = (x*cos(cp[2]) -
                            y*sin(cp[2]))
                    yp = (x*sin(cp[2]) +
                            y*cos(cp[2]))

                    #Store adhesions
                    Adh0[ad] = np.array([xp, yp]) + cell[:2]
                    Adh[ad] = np.array([xp, yp])
                    
            #sanity check to make sure all the adhesions are inside

    else:
        # When there are no adhesions, the cell
        # protrudes so that overlap energy is minimized
        theta = cp[2]
        xc, yc = cell[:2]
        aval, bval = a, b

        #create ellipses
        ell_i_in = create_ellipse((xc, yc), 
                    (aval/2, bval/2), theta)
        ell_i_out = create_ellipse((xc, yc), 
                    (aval, bval), theta)

        # Check if cell overlaps inner ellipses     
        if one_cell_overlap(ell_i_in, ell_i_out, 
                    num, cells, cparams):

            #Find the correct value of 'a'
            a_max = a
            a_min = cp[0]

            #Iterate till a value is found
            while (a_max - a_min > 1e-4):
                aval = (a_max - a_min)/2
                xc, yc = cell[:2]

                #create ellipses
                ell_i_in = create_ellipse((xc, yc), 
                    (aval/2, bval/2), theta)
                ell_i_out = create_ellipse((xc, yc), 
                    (aval, bval), theta)

                #check if this cell overlaps
                if one_cell_overlap(ell_i_in, ell_i_out, 
                    num, cells, cparams):
                    a_max = aval 
                    a_min = a_min
                else:
                    a_max = a_max
                    a_min = aval    

            #check if the cell remains in the contracted phase
            if (abs(aval - cp[0]) > 1e-4):
                #contraction phase
                cp[3] = 0
            elif cp[3] != -1:
                #no contraction phase
                cp[3] = -1
            else:
                #Don't update 'a'                        
                return cell, cp, Adh, Adh0        

            #Update a
            cp[0] = aval

            if (cparams[4*num+3] == 0):
                Adh = np.zeros((Nadh, 2)) - 1e8
                Adh0 = np.zeros((Nadh, 2)) - 1e8

                #Front and back points must be connected
                #front
                x, y = aval, 0.0
                xp = (x*cos(cp[2]) - y*sin(cp[2]))
                yp = (x*sin(cp[2]) + y*cos(cp[2]))
                Adh0[0] = np.array([xp, yp]) + cell[:2]
                Adh[0] = np.array([xp, yp])

                #back
                x, y = -aval, 0.0
                xp = (x*cos(cp[2]) - y*sin(cp[2]))
                yp = (x*sin(cp[2]) + y*cos(cp[2]))
                Adh0[1] = np.array([xp, yp]) + cell[:2]
                Adh[1] = np.array([xp, yp])

                n = 2
                for i in range(Nadh-2):
                
                    phi = rng.uniform(0, float(2*pi))
                    rho = rng.uniform(0, 1)

                    #In the system frame of reference
                    x = cp[0]*sqrt(rho)*cos(phi)
                    y = b*sqrt(rho)*sin(phi)

                    if (rng.random() < k_plus*dt):
                        #Transform into ellipse
                        xp = (x*cos(cp[2]) -
                                y*sin(cp[2]))
                        yp = (x*sin(cp[2]) +
                                y*cos(cp[2]))

                        #Store adhesions
                        Adh0[i] = np.array([xp, yp]) + cell[:+2]
                        Adh[i] = np.array([xp, yp])   

        else:
            #No adhesions and no overlaps
            #Update phase and a
            cp[0] = a
            cp[3] = 0

            #Completely new adhesions only if phase time > 3
            if (cparams[4*num+3] > 3):
                Adh = np.zeros((Nadh, 2)) - 1e8
                Adh0 = np.zeros((Nadh, 2)) - 1e8
    
                #Front and back points must be connected
                #front
                x, y = cp[0], 0.0
                xp = (x*cos(cp[2]) - y*sin(cp[2]))
                yp = (x*sin(cp[2]) + y*cos(cp[2]))
                Adh0[0] = np.array([xp, yp]) + cell[:2]
                Adh[0] = np.array([xp, yp])
    
                #back
                x, y = -cp[0], 0.0
                xp = (x*cos(cp[2]) - y*sin(cp[2]))
                yp = (x*sin(cp[2]) + y*cos(cp[2]))
                Adh0[1] = np.array([xp, yp]) + cell[:2]
                Adh[1] = np.array([xp, yp])
    
                n = 2
                for i in range(Nadh-2):
                
                    phi = rng.uniform(0, float(2*pi))
                    rho = rng.uniform(0, 1)
    
                    #In the system frame of reference
                    x = cp[0]*sqrt(rho)*cos(phi)
                    y = b*sqrt(rho)*sin(phi)
    
                    if (rng.random() < k_plus*dt):
                        #Transform into ellipse
                        xp = (x*cos(cp[2]) -
                                y*sin(cp[2]))
                        yp = (x*sin(cp[2]) +
                                y*cos(cp[2]))
    
                        #Store adhesions
                        Adh0[i] = np.array([xp, yp]) + cell[:+2]
                        Adh[i] = np.array([xp, yp])
                             
    return cell, cp, Adh, Adh0

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
