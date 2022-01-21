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
    ells_in = []
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

    a = 6
    b = 3
    L = 40
    Num = 9
    Nadh = 64
    k_plus = 0.5
    dt = 1
    lamda = 0.6
    tau = 30
    T_S = 2000
    k_s = 10.0
    k_m = 1
    k_out_out = 5
    k_in_out = k_out_out*100
    k_in_in = k_in_out*100
    fThreshold = 10.0
    k_b = 0.005
    k_f = 0.0005
    alpha = 10
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
    one_energy_vals = []
    total_energy_vals = []
    for i in np.arange(2000):
        cells[12] = x - 12 + i*eps 
        Ovlaps = Cell_funcs.find_overlaps(cells, cparams, Ovlaps)
        one_energy_vals.append(energy.overlap_energy(cells, cparams, 4, Ovlaps, k_out_out, k_in_out, k_in_in))
        total_energy_vals.append(energy.total_overlap_energy(cells, cparams, Ovlaps, k_out_out, k_in_out, k_in_in))
    cells[12] = x0   

    plt.plot(x - 12 + eps*np.arange(2000), one_energy_vals)  
    plt.title("Change in Overlap Energy of the center cell as a function of x (Area^2)")
    plt.savefig("check/one_energy_center_x_area_sqr.png")
    plt.close()  

    plt.plot(x - 12 + eps*np.arange(2000), total_energy_vals)  
    plt.title("Change in Total Overlap Energy of the center cell as a function of x (Area^2)")
    plt.savefig("check/total_energy_center_x_area_sqr.png")
    plt.close()  

    plt.plot((x - 12 + eps*np.arange(2000))[850:1150], one_energy_vals[850:1150])  
    plt.title("Change in Overlap Energy of the center cell as a function of x (Area^2 - ZOOM)")
    plt.savefig("check/one_energy_center_x_area_sqr_zoom.png")
    plt.close() 

    plt.plot((x - 12 + eps*np.arange(2000))[850:1150], total_energy_vals[850:1150])    
    plt.title("Change in Total Overlap Energy of the center cell as a function of x (Area^2 - ZOOM)")
    plt.savefig("check/total_energy_center_x_area_sqr_zoom.png")
    plt.close()

    eps = 15/2000
    y0 = cells[13] 
    y = cells[13]
    one_energy_vals = []
    total_energy_vals = []
    for i in np.arange(2000):
        cells[13] = y - 7.5 + i*eps 
        Ovlaps = Cell_funcs.find_overlaps(cells, cparams, Ovlaps)
        one_energy_vals.append(energy.overlap_energy(cells, cparams, 4, Ovlaps, k_out_out, k_in_out, k_in_in))
        total_energy_vals.append(energy.total_overlap_energy(cells, cparams, Ovlaps, k_out_out, k_in_out, k_in_in))
    cells[13] = y0

    plt.plot(y - 7.5 + eps*np.arange(2000), one_energy_vals)  
    plt.title("Overlap Energy vs y of the center cell (Area^2)")
    plt.savefig("check/one_energy_center_y_area_sqr.png")
    plt.close()  

    plt.plot(y - 7.5 + eps*np.arange(2000), total_energy_vals)  
    plt.title("Total Overlap Energy vs y of the center cell (Area^2)")
    plt.savefig("check/total_energy_center_y_area_sqr.png")
    plt.close()  

    plt.plot((y - 7.5 + eps*np.arange(2000))[850:1150], one_energy_vals[850:1150])  
    plt.title("Overlap Energy vs y of the center cell (Area^2 - ZOOM)")
    plt.savefig("check/one_energy_center_y_area_sqr_zoom.png")
    plt.close() 

    plt.plot((y - 7.5 + eps*np.arange(2000))[850:1150], total_energy_vals[850:1150])    
    plt.title("Total Overlap Energy vs y of the center cell (Area^2 - ZOOM)")
    plt.savefig("check/total_energy_center_y_area_sqr_zoom.png")
    plt.close()

    eps = (2*np.pi)/2000
    th0 = cells[14]
    th = cells[14]
    one_energy_vals = []
    total_energy_vals = []
    for i in np.arange(2000):
        cells[14] = th - np.pi + i*eps 
        Ovlaps = Cell_funcs.find_overlaps(cells, cparams, Ovlaps)
        one_energy_vals.append(energy.overlap_energy(cells, cparams, 4, Ovlaps, k_out_out, k_in_out, k_in_in))
        total_energy_vals.append(energy.total_overlap_energy(cells, cparams, Ovlaps, k_out_out, k_in_out, k_in_in))
    cells[14] = th0

    plt.plot(th - np.pi + eps*np.arange(2000), one_energy_vals)  
    plt.title("Overlap Energy vs theta of the center cell (Area^2)")
    plt.savefig("check/one_energy_center_theta_area_sqr.png")
    plt.close()  

    plt.plot(th - np.pi + eps*np.arange(2000), total_energy_vals)  
    plt.title("Total Overlap Energy vs theta of the center cell (Area^2)")
    plt.savefig("check/total_energy_center_theta_area_sqr.png")
    plt.close()

    plt.plot(th - np.pi + eps*np.arange(2000), one_energy_vals)  
    plt.title("Overlap Energy vs theta of the center cell (Area^2 - ZOOM)")
    plt.savefig("check/one_energy_center_theta_area_sqr_zoom.png")
    plt.close()  

    plt.plot(th - np.pi + eps*np.arange(2000), total_energy_vals)  
    plt.title("Total Overlap Energy vs theta of the center cell (Area^2 - ZOOM)")
    plt.savefig("check/total_energy_center_theta_area_sqr_zoom.png")
    plt.close()