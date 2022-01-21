from matplotlib import pyplot as plt
from shapely.geometry import Point, LineString
from shapely import affinity
from shapely.ops import unary_union, polygonize
from itertools import combinations
from matplotlib.patches import Polygon
from matplotlib import collections  as mc
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

def plot_ellipses(cells, cparams, Adh, Adh0, t):

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
                        edgecolor="black", linewidth = 2)
        ax.add_patch(patch)

    #ells_in = unary_union(ells_in)
    for e in polygonize(ells_in):
        verts = np.array(e.exterior.coords.xy)
        patch = Polygon(verts.T, color = 'red', alpha = 0.8)
        ax.add_patch(patch)

    #Polarity vector
    ax.quiver(*np.array([X, Y]), np.array(U), np.array(V), units='xy', angles="xy", scale=1)

    lc = mc.LineCollection(liness, colors="black", linewidths=1)
    ax.add_collection(lc)
    x1 = [j[0][0] for j in liness]
    y1 = [j[0][1] for j in liness]
    x2 = [j[1][0] for j in liness]
    y2 = [j[1][1] for j in liness]
    ax.scatter(x1, y1, color="orange", s=4)
    ax.scatter(x2, y2, color="green", s=4)

    multiplier = 5.0
    for axis, setter in [(ax.xaxis, ax.set_xlim), (ax.yaxis, ax.set_ylim)]:
        vmin, vmax = axis.get_data_interval()
        vmin = multiplier * np.floor(vmin / multiplier)
        vmax = multiplier * np.ceil(vmax / multiplier)
        setter([vmin, vmax])

    plt.title("x = {:.4f}, y = {:.4f}, theta = {:.4f}".format(cells[0], cells[1], cparams[2]))
    plt.savefig("plots/{}.png".format(t))
    plt.close()

def plot_from_data():

    DATA = np.load("DATA.npz")
    CELLS = DATA["CELLS"]
    CPARAMS = DATA["CPARAMS"]
    ADH = DATA["ADH"]
    ADH0 = DATA["ADH0"]

    for t in range(CELLS.shape[0]):

        plot_ellipses(CELLS[t], CPARAMS[t], 
                ADH[t], ADH0[t], t)

if __name__ == "__main__":
    plot_from_data()