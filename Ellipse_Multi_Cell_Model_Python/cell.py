from shapely.geometry.point import Point
from shapely import affinity
import numpy as np
from math import sin, cos, pi, exp, sqrt

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
        ell2 = create_ellipse((cells[3*ind:(3*ind+2)]), (a,b),
                               cells[3*ind+2])

        if ell1.boundary.intersects(ell2.boundary):
            return False

    return True

def random_cells(L, a, b, Num, rng):

    # x y dtheta
    cells = np.zeros(Num*3)
    # a b theta phase
    cparams = np.zeros(Num*4)
    Ovlaps = np.zeros((Num, Num)) - 1e8

    ind = 0
    while (ind < Num):

        x, y, theta = ([rng.uniform(0, L), rng.uniform(0, L),
                        rng.uniform(0, 360)])

        # check for collisions
        if noCollision(cells, ind, a, b, x, y, theta):
            cells[3*ind:(3*ind+3)] = np.array([x, y, 0])
            cparams[4*ind:(4*ind+4)] = np.array([a, b, theta, 0])
            ind += 1

    return cells, cparams, Ovlaps

def find_overlaps(cells, cparams, Ovlaps):

    for i in range(cells.shape[0]//3):
        ell_i = create_ellipse((cells[3*i],cells[3*i+1]), (
                                cparams[4*i],cparams[4*i+1]), cparams[4*i+2])
        for j in range(i+1, cells.shape[0]//3):
            ell_j = create_ellipse((cells[3*j],cells[3*j+1]), (
                                    cparams[4*j],cparams[4*j+1]), cparams[4*j+2])

            if ell_i.boundary.intersects(ell_j.boundary):
                Ovlaps[i][j] = 1
                Ovlaps[j][i] = 1

    return Ovlaps
