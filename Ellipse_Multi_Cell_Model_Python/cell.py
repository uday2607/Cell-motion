from shapely.geometry.point import Point
from shapely import affinity
import numpy as np
from math import sin, cos, pi, exp

def create_ellipse(center, lengths, angle=0):
    """
    create a shapely ellipse. adapted from
    https://gis.stackexchange.com/a/243462
    """
    circ = Point(center).buffer(1, resolution=1000)
    ell = affinity.scale(circ, int(lengths[0]), int(lengths[1]))
    ellr = affinity.rotate(ell, angle)
    return ellr

def noCollision(cells, ind, a, b, x, y, theta):

    if ind == 0:
        return True

    ell1 = create_ellipse((x,y), (a,b), theta)

    for i in range(ind):
        ell2 = create_ellipse((cell[3*ind:3*ind+2]), (a,b), cell[3*ind+2])

        if ell1.boundary.intersects(ell2.boundary):
            return False

    return True

def random_cells(L, a, b, Num, rng):

    # x y theta
    cells = np.zeros(Num*3)
    # a b phase
    cparams = np.zeros(Num*3)
    Ovlaps = np.zeros((Num, Num)) - 1

    ind = 0
    while (ind < Num):

        x, y, theta = ([rng.uniform(0, L), rng.uniform(0, L),
                        rng.uniform(0, 360)])

        # check for collisions
        if noCollision(cells, ind, a, b, x, y, theta):
            cells[3*ind:3*ind+3] = np.array([x, y, theta])
            cparams[3*ind:3*ind+3] = np.array([a, b, 0])
            ind += 1

    return cells, cparams, Ovlaps
