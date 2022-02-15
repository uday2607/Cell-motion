import numpy as np
import numba as nb
from math import sin, cos, tan, pi, exp, sqrt

"""Polygon Area function"""
#From: https://stackoverflow.com/a/58515054
@nb.jit(nopython = True, nogil = True)
def polygon_area(x_y):
    x = x_y[:,0]
    y = x_y[:,1]

    S1 = np.sum(x*np.roll(y,-1))
    S2 = np.sum(y*np.roll(x,-1))

    area = .5*np.absolute(S1 - S2)

    return area
""""""    