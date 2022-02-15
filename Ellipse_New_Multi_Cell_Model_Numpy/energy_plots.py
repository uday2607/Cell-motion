import Cell_funcs
import Adh_funcs
import energy
import numpy as np
import math

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

def plot_ellipses(cells, cparams, title):

    ells_out = []
    ells_mid = []
    ells_in = []
    liness = []
    X = []
    Y = []
    U = []
    V = []
    C = [1]*(cells.shape[0]//3)

    for i in range(cells.shape[0]//3):
        temp_ell = create_ellipse(cells[3*i:3*i+2],
                    (cparams[4*i]*4,cparams[4*i+1]*4),cparams[4*i+2])
        ells_out.append(temp_ell)
        temp_ell = create_ellipse(cells[3*i:3*i+2],
                    (cparams[4*i]*2,cparams[4*i+1]*2),cparams[4*i+2])
        ells_mid.append(temp_ell)
        temp_ell = create_ellipse(cells[3*i:3*i+2],
                    (cparams[4*i],cparams[4*i+1]),cparams[4*i+2])
        ells_in.append(temp_ell)

        #Polarity vector
        X.append(cells[3*i])
        Y.append(cells[3*i+1])
        U.append(cparams[4*i]*math.cos(cparams[4*i+2]))
        V.append(cparams[4*i]*math.sin(cparams[4*i+2]))

    fig, ax = plt.subplots(subplot_kw={'aspect': 'equal'})
    for e in polygonize(ells_out):
        verts = np.array(e.exterior.coords.xy)
        patch = Polygon(verts.T, facecolor = 'violet', alpha = 0.2, linestyle = "-",
                        edgecolor="black", linewidth = 2, zorder=30)
        ax.add_patch(patch)

    for e in polygonize(ells_mid):
        verts = np.array(e.exterior.coords.xy)
        patch = Polygon(verts.T, facecolor = 'blue', alpha = 0.5, linestyle = "-",
                        edgecolor="black", linewidth = 1.5, zorder=30)
        ax.add_patch(patch)    

    for e in polygonize(ells_in):
        verts = np.array(e.exterior.coords.xy)
        patch = Polygon(verts.T, color = 'red', alpha = 0.8, zorder=50)
        ax.add_patch(patch)

    #Polarity vector
    ax.quiver(*np.array([X, Y]), np.array(U), np.array(V), units='xy', angles="xy", scale=1)

    multiplier = 10.0
    for axis, setter in [(ax.xaxis, ax.set_xlim), (ax.yaxis, ax.set_ylim)]:
        vmin, vmax = axis.get_data_interval()
        vmin = multiplier * np.floor(vmin / multiplier)
        vmax = multiplier * np.ceil(vmax / multiplier)
        setter([vmin, vmax])

    plt.title("x = {:.4f}, y = {:.4f}, theta = {:.4f}".format(cells[0], cells[1], cparams[2]))
    plt.savefig("check/"+title+".png")
    plt.close()

if __name__ == '__main__':

    a = 4
    b = 2
    L = 1000
    Num = 9
    Nadh = 64
    k_plus = 0.025
    lamda = 0.2
    tau = 30
    dt = 60/tau
    T_S_c = 1000
    T_S_p = 1000
    k_s = 40.0
    k_m = 1
    E_const = np.zeros((3, 3))
    E_const[2, 2] = -0.0025
    E_const[1, 1] = 1.0
    E_const[1, 0] = 5000.0
    E_const[0, 1] = 5000.0
    E_const[0, 0] = 10000.0
    fThreshold = 40.0/k_s
    k_b = 0.05
    k_f = 0.005
    alpha = 25
    a_min = a
    for i in range(tau):
        a_min -= a_min*lamda*dt/tau
    rng = np.random.default_rng()

    # Spawn Cells and Adhesions
    #cells, cparams, Ovlaps = Cell_funcs.random_cells(L, a, b, Num)
    #cells, cparams, Ovlaps = Cell_funcs.preset(a, b, Num)
    cells, cparams, Ovlaps = Cell_funcs.custom_preset(a, b, Num)
    cparams0 = cparams.copy()

    #Find overlap between cells
    Ovlaps = Cell_funcs.find_overlaps(cells, cparams, Ovlaps)

    plot_ellipses(cells, cparams, "system")

    # For middle cell
    eps = 24/2000
    x0 = cells[12]
    x = cells[12]
    total_energy_vals = []
    for i in np.arange(2000):
        cells[12] = x - 12 + i*eps 
        Ovlaps = Cell_funcs.find_overlaps(cells, cparams, Ovlaps)
        total_energy_vals.append(energy.total_overlap_energy(cells, cparams, Ovlaps, E_const))
    cells[12] = x0

    plt.plot(x - 12 + eps*np.arange(2000), total_energy_vals)  
    plt.title("Total Overlap Energy vs x of the center cell (Area^2)")
    plt.savefig("check/total_energy_center_x_area_sqr.png")
    plt.close()

    