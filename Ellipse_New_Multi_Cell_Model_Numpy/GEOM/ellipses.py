import numpy as np
import numba as nb
from math import sin, cos, tan, pi, exp, sqrt
from GEOM.intersection import *
from GEOM.area import *

#Function to create an ellipse using pure numpy and numba
@nb.jit(nopython = True, nogil = True)
def create_ellipse(center, lengths, angle=0):

    #Create empty array
    N = 45
    phi = -np.arange(N)*2*pi/N
    ells = np.empty((3, N, 2))

    #Cell boundary coordinates
    x = lengths[0]*np.cos(phi)
    y = lengths[1]*np.sin(phi)

    #Transform the ellipse
    ells[0, :, 0] = x*np.cos(angle) - y*np.sin(angle) + center[0]
    ells[0, :, 1] = x*np.sin(angle) + y*np.cos(angle) + center[1]

    #Transform the ellipse
    ells[1, :, 0] = 2.0*x*np.cos(angle) - 2.0*y*np.sin(angle) + center[0]
    ells[1, :, 1] = 2.0*x*np.sin(angle) + 2.0*y*np.cos(angle) + center[1]

    #Transform the ellipse
    ells[2, :, 0] = 4.0*x*np.cos(angle) - 4.0*y*np.sin(angle) + center[0]
    ells[2, :, 1] = 4.0*x*np.sin(angle) + 4.0*y*np.cos(angle) + center[1]

    return ells

# check if a point is in the ellipse
@nb.jit(nopython = True, nogil = True)
def in_ellipse(x, y, xc, yc, a, b, theta):

    if (((x-xc)*cos(theta)+(y-yc)*sin(theta))**2/a**2 +
        ((x-xc)*sin(theta)-(y-yc)*cos(theta))**2/b**2) < 1:
        return 1
    else:
        return 0    

# Check if two ellipses intersect
@nb.jit(nopython = True, nogil = True)
def check_ell_intersection(x1, y1, a1, b1, t1,
                         x2, y2, a2, b2, t2):

    inside = 0
    for i in range(45):
        phi = 2*pi*(44-i)/45
        x_ = a1*cos(phi)
        y_ = b1*sin(phi)

        # coordinate of the first ellipse
        x = x_*cos(t1) - y_*sin(t1) + x1 
        y = x_*sin(t1) + y_*cos(t1) + y1

        # check if the point is in ellipse 2
        xt = (x - x2)*cos(t2) + (y - y2)*sin(t2)
        yt = -(x - x2)*sin(t2) + (y - y2)*cos(t2)
        if xt**2/a2**2 + yt**2/b2**2 < 1:
            inside = 1
            break        

    return inside

# Calculates the overlap area of the ellipses
@nb.jit(nopython = True, nogil = True)
def ell_ell_area(ell1, ell2):
    
    intersection = clip(ell1, ell2)
    if intersection.shape[0] != 0:
        area = polygon_area(intersection)
    else:
        area = 0.0    

    return area