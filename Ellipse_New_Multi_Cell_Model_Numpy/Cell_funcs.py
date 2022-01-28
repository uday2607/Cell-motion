from ellipses import *
from random_funcs import *

@nb.jit(nopython = True, nogil = True)
def noCollision(cells, cparams, ind, a, b, x, y, theta):

    if ind == 0:
        return True

    x1, y1 = x, y
    a1, b1 = a, b
    t1 = theta
    for i in range(ind):
        x2, y2 = cells[3*i:(3*i+2)]
        a2, b2 = a, b
        t2 = cparams[4*i+2]

        if check_ell_intersection(x1, y1, a1, b1, t1,
                x2, y2, a2, b2, t2):
            return 0

    return 1

@nb.jit(nopython = True, nogil = True)
def preset(a, b, Num):

   # x y dtheta
    cells = np.zeros(Num*3)
    # a b theta phase
    cparams = np.zeros(Num*4)
    Ovlaps = np.zeros((Num, Num)) - 1e8

    for ind in range(Num):

        if ind // (Num//2) == 0 & ind != Num//2:
            theta = pi
        else:
            theta = 0

        cells[3*ind:(3*ind+3)] = np.array([20+2*(a-a/4)*ind, 10, 0])
        cparams[4*ind:(4*ind+4)] = np.array([a, b, theta, 1])

    return cells, cparams, Ovlaps

@nb.jit(nopython = True, nogil = True)
def linear_preset(a, b, Num):

    # x y dtheta
    cells = np.zeros(Num*3)
    # a b theta phase
    cparams = np.zeros(Num*4)
    Ovlaps = np.zeros((Num, Num)) - 1e8

    for ind in range(Num):
        cells[3*ind:(3*ind+3)] = np.array([20+(a-a/4)*ind, 10, 0])
        cparams[4*ind:(4*ind+4)] = np.array([a/2, b, 0, 1])

    return cells, cparams, Ovlaps

@nb.jit(nopython = True, nogil = True)
def two_linear_preset(a, b, Num):

    num = Num//2

    # x y dtheta
    cells = np.zeros(Num*3)
    # a b theta phase
    cparams = np.zeros(Num*4)
    Ovlaps = np.zeros((Num, Num)) - 1e8

    for ind in range(Num):
        cells[3*ind:(3*ind+3)] = np.array([20+2*(a-a/4)*(ind%num), 10+2*(b-b/4)*(ind//num), 0])
        cparams[4*ind:(4*ind+4)] = np.array([a, b, 0, 1])

    return cells, cparams, Ovlaps

@nb.jit(nopython = True, nogil = True)
def custom_preset(a, b, Num):

    Num = int(sqrt(Num))**2
    num = int(sqrt(Num))

    # x y dtheta
    cells = np.zeros(Num*3)
    # a b theta phase
    cparams = np.zeros(Num*4)
    Ovlaps = np.zeros((Num, Num)) - 1e8

    for ind in range(Num):
        cells[3*ind:(3*ind+3)] = np.array([20+2*(a-1)*(ind%num), 10+2*(b-0.5)*(ind//num), 0])
        cparams[4*ind:(4*ind+4)] = np.array([a, b, 0, 1])

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
        if noCollision(cells, cparams, ind, a, b, x, y, theta):
            cells[3*ind:(3*ind+3)] = np.array((x, y, 0))
            cparams[4*ind:(4*ind+4)] = np.array((a, b, theta, 1))
            ind += 1

    return cells, cparams, Ovlaps

@nb.jit(nopython = True, nogil = True)
def find_overlaps(cells, cparams, Ovlaps):

    for i in range(cells.shape[0]//3):
        x1, y1 = cells[3*i], cells[3*i+1]
        a1, b1 = cparams[4*i], cparams[4*i+1]
        t1 = cparams[4*i+2]
        for j in range(i+1, cells.shape[0]//3):
            x2, y2 = cells[3*j], cells[3*j+1]
            a2, b2 = cparams[4*j], cparams[4*j+1]
            t2 = cparams[4*j+2]

            if check_ell_intersection(x1, y1, a1, b1, t1,
                                    x2, y2, a2, b2, t2):
                Ovlaps[i][j] = 1
                Ovlaps[j][i] = 1
            else:
                Ovlaps[i][j] = 0
                Ovlaps[j][i] = 0

    return Ovlaps

@nb.jit(nopython = True, nogil = True)
def find_overlaps_ind(cells, cparams, Ovlaps, i):

    x1, y1 = cells[3*i], cells[3*i+1]
    a1, b1 = cparams[4*i], cparams[4*i+1]
    t1 = cparams[4*i+2]
    for j in range(i+1, cells.shape[0]//3):
        x2, y2 = cells[3*j], cells[3*j+1]
        a2, b2 = cparams[4*j], cparams[4*j+1]
        t2 = cparams[4*j+2]

        if check_ell_intersection(x1, y1, a1, b1, t1,
                                x2, y2, a2, b2, t2):
            Ovlaps[i][j] = 1
            Ovlaps[j][i] = 1
        else:
            Ovlaps[i][j] = 0
            Ovlaps[j][i] = 0

    return Ovlaps
