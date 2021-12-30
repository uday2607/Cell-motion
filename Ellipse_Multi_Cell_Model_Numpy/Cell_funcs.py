from ellipses import *
from random_funcs import *

@nb.jit(nopython = True, nogil = True)
def noCollision(cells, ind, a, b, x, y, theta):

    if ind == 0:
        return True

    ell1 = create_ellipse(np.array((x,y)), np.array((a,b)), theta)

    for i in range(ind):
        ell2 = create_ellipse(cells[3*ind:(3*ind+2)], np.array((a,b)),
                               cells[3*ind+2])

        intersection = ellipse_intersection(ell1, ell2)
        if intersection.shape[0] == 0:
            return 0

    return 1

@nb.jit(nopython = True, nogil = True)
def preset(a, b, Num):
    # x y dtheta
    cells = np.zeros(Num*3)
    # a b theta phase
    cparams = np.zeros(Num*4)
    Ovlaps = np.zeros((Num, Num)) - 1e8

    ind = 1
    x = 0
    y = 0
    theta = 0
    cells[3*ind:(3*ind+3)] = np.array([x, y, 0])
    cparams[4*ind:(4*ind+4)] = np.array([a, b, theta, 1])

    ind = 0
    x = -2
    y = -7
    theta = 45*pi/180
    cells[3*ind:(3*ind+3)] = np.array([x, y, 0])
    cparams[4*ind:(4*ind+4)] = np.array([a, b, theta, 1])

    return cells, cparams, Ovlaps

@nb.jit(nopython = True, nogil = True)
def random_cells(L, a, b, Num):

    # x y dtheta
    cells = np.zeros(Num*3)
    # a b theta phase
    cparams = np.zeros(Num*4)
    Ovlaps = np.zeros((Num, Num)) - 1e8

    ind = 0
    while (ind < Num):
        x = uniform_double(0, L, 1)[0]
        y = uniform_double(0, L, 1)[0]
        theta = uniform_double(0, 2*pi, 1)[0]

        # check for collisions
        if noCollision(cells, ind, a, b, x, y, theta):
            cells[3*ind:(3*ind+3)] = np.array((x, y, 0))
            cparams[4*ind:(4*ind+4)] = np.array((a, b, theta, 1))
            ind += 1

    return cells, cparams, Ovlaps

@nb.jit(nopython = True, nogil = True)
def find_overlaps(cells, cparams, Ovlaps):

    for i in range(cells.shape[0]//3):
        ell_i = create_ellipse((cells[3*i],cells[3*i+1]), (
                                cparams[4*i],cparams[4*i+1]), cparams[4*i+2])
        for j in range(i+1, cells.shape[0]//3):
            ell_j = create_ellipse((cells[3*j],cells[3*j+1]), (
                                    cparams[4*j],cparams[4*j+1]), cparams[4*j+2])

            if ellipse_intersection(ell_i, ell_j).shape[0] != 0:
                Ovlaps[i][j] = 1
                Ovlaps[j][i] = 1

    return Ovlaps
