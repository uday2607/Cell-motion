from matplotlib import pyplot as plt
from shapely.geometry import Point, LineString
from shapely import affinity
from shapely.ops import unary_union, polygonize
from itertools import combinations
from matplotlib.patches import Polygon
import numpy as np
import math

def create_ellipse(center, lengths, angle=0):
    """
    create a shapely ellipse. adapted from
    https://gis.stackexchange.com/a/243462
    """
    circ = Point(center).buffer(1, resolution=1000)
    ell = affinity.scale(circ, int(lengths[0]), int(lengths[1]))
    ellr = affinity.rotate(ell, angle)
    return ellr

def plot_ellipses(cells, cparams, a, b, t):

    ells_out = []
    ells_in = []
    x_min = 0
    x_max = 0
    y_min = 0
    y_max = 0

    for i in range(cells.shape[0]//3):
        temp_ell = create_ellipse(cells[3*i:3*i+2],
                    (cparams[4*i],cparams[4*i+1]),cparams[4*i+2])
        ells_out.append(temp_ell)
        temp_ell = create_ellipse(cells[3*i:3*i+2],
                    (cparams[4*i]/2,cparams[4*i+1]/2),cparams[4*i+2])
        ells_in.append(temp_ell)

        if (x_min > cells[3*i]):
            x_min = cells[3*i]
        if (x_max <= cells[3*i]):
            x_max = cells[3*i]
        if (y_min > cells[3*i+1]):
            y_min = cells[3*i+1]
        if (y_max <= cells[3*i+1]):
            y_max = cells[3*i+1]


    fig, ax = plt.subplots(subplot_kw={'aspect': 'equal'})
    ells_out = unary_union(ells_out)
    for e in polygonize(ells_out):
        verts = np.array(e.exterior.coords.xy)
        patch = Polygon(verts.T, color = 'blue', alpha = 0.5)
        ax.add_patch(patch)

    ells_in = unary_union(ells_in)
    for e in polygonize(ells_in):
        verts = np.array(e.exterior.coords.xy)
        patch = Polygon(verts.T, color = 'red', alpha = 0.8)
        ax.add_patch(patch)

    ax.set_xlim(-2*a+x_min, x_max+2*a)
    ax.set_ylim(-2*b+y_min, y_max+2*b)

    plt.savefig("plots/{}.png".format(t))
    plt.close()
