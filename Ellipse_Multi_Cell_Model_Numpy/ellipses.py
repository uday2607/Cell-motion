from matplotlib.pyplot import boxplot
import numpy as np
import numba as nb
from math import sin, cos, tan, pi, exp, sqrt
from sys import exit

#Function to create an ellipse using pure numpy and numba
@nb.jit(nopython = True, nogil = True)
def create_ellipse(center, lengths, angle=0):

    #Create empty array
    N = 45
    phi = -np.arange(N)*2*pi/N
    ell_out = np.empty((N, 2))
    ell_in = np.empty((N, 2))

    #Cell boundary coordinates
    x = lengths[0]*np.cos(phi)
    y = lengths[1]*np.sin(phi)

    #Transform the ellipse
    ell_out[:, 0] = x*np.cos(angle) - y*np.sin(angle) + center[0]
    ell_out[:, 1] = x*np.sin(angle) + y*np.cos(angle) + center[1]

    #Transform the ellipse
    ell_in[:, 0] = 0.5*x*np.cos(angle) - 0.5*y*np.sin(angle) + center[0]
    ell_in[:, 1] = 0.5*x*np.sin(angle) + 0.5*y*np.cos(angle) + center[1]

    return ell_out, ell_in

"""Polygon Clipper functions"""
#From: https://github.com/mdabdk/sutherland-hodgman
#Helper functions
@nb.jit(nopython = True, nogil = True)
def is_inside(p1, p2, q):
    R = (p2[0] - p1[0]) * (q[1] - p1[1]) - (p2[1] - p1[1]) * (q[0] - p1[0])
    if R <= 0:
        return 1
    else:
        return 0

@nb.jit(nopython = True, nogil = True)
def compute_intersection(p1, p2, p3, p4):

    # if first line is vertical
    if p2[0] - p1[0] == 0:
        x = p1[0]
        
        # slope and intercept of second line
        m2 = (p4[1] - p3[1]) / (p4[0] - p3[0])
        b2 = p3[1] - m2 * p3[0]
        
        # y-coordinate of intersection
        y = m2 * x + b2
    
    # if second line is vertical
    elif p4[0] - p3[0] == 0:
        x = p3[0]
        
        # slope and intercept of first line
        m1 = (p2[1] - p1[1]) / (p2[0] - p1[0])
        b1 = p1[1] - m1 * p1[0]
        
        # y-coordinate of intersection
        y = m1 * x + b1
    
    # if neither line is vertical
    else:
        m1 = (p2[1] - p1[1]) / (p2[0] - p1[0])
        b1 = p1[1] - m1 * p1[0]
        
        # slope and intercept of second line
        m2 = (p4[1] - p3[1]) / (p4[0] - p3[0])
        b2 = p3[1] - m2 * p3[0]
    
        # x-coordinate of intersection
        x = (b2 - b1) / (m1 - m2)
    
        # y-coordinate of intersection
        y = m1 * x + b1
    
    intersection = np.empty((1, 2))
    intersection[0, 0] = x 
    intersection[0, 1] = y
    
    return intersection

#Clipper function
@nb.jit(nopython = True, nogil = True)
def clip(subject_polygon,clipping_polygon):
        
    final_polygon = subject_polygon.copy()
    
    for i in range(clipping_polygon.shape[0]):
        
        # stores the vertices of the next iteration of the clipping procedure
        next_polygon = final_polygon.copy()
        
        # stores the vertices of the final clipped polygon
        final_polygon = np.empty((0, 2))
        
        # these two vertices define a line segment (edge) in the clipping
        # polygon. It is assumed that indices wrap around, such that if
        # i = 1, then i - 1 = K.
        c_edge_start = clipping_polygon[i - 1]
        c_edge_end = clipping_polygon[i]
        
        for j in range(next_polygon.shape[0]):
            
            # these two vertices define a line segment (edge) in the subject
            # polygon
            s_edge_start = next_polygon[j - 1]
            s_edge_end = next_polygon[j]
            
            if is_inside(c_edge_start,c_edge_end,s_edge_end):
                if not is_inside(c_edge_start,c_edge_end,s_edge_start):
                    intersection = compute_intersection(s_edge_start,s_edge_end,c_edge_start,c_edge_end)
                    final_polygon = np.vstack((final_polygon, intersection))
                s_edge_end_ = np.empty((1, 2))
                s_edge_end_[0] = s_edge_end    
                final_polygon = np.vstack((final_polygon, s_edge_end_))
            elif is_inside(c_edge_start,c_edge_end,s_edge_start):
                intersection = compute_intersection(s_edge_start,s_edge_end,c_edge_start,c_edge_end)
                final_polygon = np.vstack((final_polygon, intersection))
    
    return final_polygon

@nb.jit(nopython = True, nogil = True)
def check_ell_intersection(x1, y1, a1, b1, t1,
                         x2, y2, a2, b2, t2):

    ell, _ = create_ellipse((x1, y1), (a1, b1), t1)
    X, Y = ell[:, 0], ell[:, 1]

    xt = (X - x2)*cos(t2) + (Y - y2)*sin(t2)
    yt = -(X - x2)*sin(t2) + (Y - y2)*cos(t2)
    check = xt**2/a2**2 + yt**2/b2**2
    inside = np.any(check <= 1)

    return inside
""""""

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

@nb.jit(nopython = True, nogil = True)
def ell_ell_area(ell1, ell2):
    
    intersection = clip(ell1, ell2)
    if intersection.shape[0] != 0:
        area = polygon_area(intersection)
    else:
        area = 0.0    

    return area    
""""""