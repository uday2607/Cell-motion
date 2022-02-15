from Cell_funcs import *
from scipy.optimize import fsolve

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
def protrusion(cells_, num, Adh_, Adh0_, cparams_, Ovlaps, dt, tau, a, a_min,
            E_const, k_s, T_S):

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
        
    # Update a and phase
    cparams[4*num] = a
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
        theta = cparams[4*num+2]
        args = (cparams[4*num], cparams[4*num+1], theta,
                Adh[num, ad, 0], Adh[num, ad, 1], cells[3*num],
                cells[3*num+1])
        xsol = fsolve(solve_center, cells[3*num:3*num+2]+np.ones(2), args=args)
        x_, y_ = xsol

        # distance to protrude
        cells[3*num], cells[3*num+1] = x_, y_

    #Change phase if the protrusion phase is over
    if cparams[4*num+3] < 0 and cparams[4*num] >= a:
        cparams[4*num+3] = 1 # +ve times -> Contraction

    return cells[3*num:3*num+3], cparams[4*num:(4*num+4)], p_flag