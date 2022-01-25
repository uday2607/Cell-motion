from matplotlib import pyplot as plt
from shapely.geometry import Point, LineString
from shapely import affinity
from shapely.ops import unary_union, polygonize
from itertools import combinations
from matplotlib.patches import Polygon
from matplotlib import collections  as mc
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import numpy as np
import math

def create_ellipse(center, lengths, angle=0):
    """
    create a shapely ellipse. adapted from
    https://gis.stackexchange.com/a/243462
    """
    circ = Point(center).buffer(1, resolution=1000)
    ell = affinity.scale(circ, lengths[0], lengths[1])
    ellr = affinity.rotate(ell, angle*180/math.pi)
    return ellr

def plot_ellipses(cells, cparams, Adh, Adh0, a, b, Y_MAX, Y_MIN, t):

    ells_out = []
    ells_in = []
    liness = []
    X = []
    Y = []
    U = []
    V = []
    C = [1]*(cells.shape[0]//3)

    for i in range(cells.shape[0]//3):
        temp_ell = create_ellipse(cells[3*i:3*i+2],
                    (cparams[4*i],cparams[4*i+1]),cparams[4*i+2])
        ells_out.append(temp_ell)
        temp_ell = create_ellipse(cells[3*i:3*i+2],
                    (cparams[4*i]/2,cparams[4*i+1]/2),cparams[4*i+2])
        ells_in.append(temp_ell)

        #Polarity vector
        X.append(cells[3*i])
        Y.append(cells[3*i+1])
        U.append(cparams[4*i]*math.cos(cparams[4*i+2]))
        V.append(cparams[4*i]*math.sin(cparams[4*i+2]))

        ind = np.arange(Adh.shape[1])[np.all(Adh[i]!=-1e8,axis=1)]

        for n in ind:
            x2, y2 = Adh[i, n, 0], Adh[i, n, 1]
            x2 += cells[3*i]
            y2 += cells[3*i+1]
            liness.append([(x2, y2), (Adh0[i, n, 0], Adh0[i, n, 1])])


    fig, ax = plt.subplots(subplot_kw={'aspect': 'equal'})
    #ells_out = unary_union(ells_out)
    for e in polygonize(ells_out):
        verts = np.array(e.exterior.coords.xy)
        patch = Polygon(verts.T, facecolor = 'blue', alpha = 0.5, linestyle = "-",
                        edgecolor="black", linewidth = 3, zorder=30)
        ax.add_patch(patch)

    #ells_in = unary_union(ells_in)
    for e in polygonize(ells_in):
        verts = np.array(e.exterior.coords.xy)
        patch = Polygon(verts.T, color = 'red', alpha = 0.8, zorder=50)
        ax.add_patch(patch)

    #Polarity vector
    ax.quiver(*np.array([X, Y]), np.array(U), np.array(V), units='xy', angles="xy", scale=1, linewidth=.5, zorder=100)

    lc = mc.LineCollection(liness, colors="black", linewidths=1.5, zorder=50)
    ax.add_collection(lc)
    x1 = [j[0][0] for j in liness]
    y1 = [j[0][1] for j in liness]
    x2 = [j[1][0] for j in liness]
    y2 = [j[1][1] for j in liness]
    ax.scatter(x1, y1, color="orange", s=4, zorder=50)
    ax.scatter(x2, y2, color="green", s=4, zorder=50)

    multiplier = 10.0
    for axis, setter in [(ax.xaxis, ax.set_xlim), (ax.yaxis, ax.set_ylim)]:
        vmin, vmax = axis.get_data_interval()
        vmin = multiplier * np.floor(vmin / multiplier)
        vmax = multiplier * np.ceil(vmax / multiplier)
        setter([vmin, vmax])
    plt.axhline(y=Y_MAX, color='r', linestyle='-')
    plt.axhline(y=Y_MIN, color='r', linestyle='-')  

    plt.tight_layout()
    #ax.xaxis.set_major_locator(MultipleLocator(1))
    #ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    #ax.xaxis.set_minor_locator(MultipleLocator(0.25))
    #ax.yaxis.set_major_locator(MultipleLocator(1))
    #ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
    #ax.yaxis.set_minor_locator(MultipleLocator(0.25))
    #plt.grid(b=True, which='both', color='black', linestyle='-', zorder=1, alpha=0.4)
    plt.title("x = {:.4f}, y = {:.4f}, theta = {:.4f}".format(cells[0], cells[1], cparams[2]))
    plt.savefig("plots/{}.png".format(t))
    plt.close()
