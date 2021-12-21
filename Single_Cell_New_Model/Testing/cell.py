import numpy as np
import numba as nb
import shapely
from shapely.geometry import Point, LineString
from matplotlib import pyplot as plt
from matplotlib.patches import Polygon
from math import cos, sin, pi

def incell(x, y, a, b, d):

    if (((x-d)/a)**2 + (y/b)**2*(
    1 + (2*d*x-d**2)/a**2)) < 1:
        return True

    return False

# params
a = 6
b = 4
d = 1
xc = 0.0
yc = 0.0
angle = 0.0
lamda = 0.4
ks = 0.1
max_num_adh = 64
k_mf = 0.02
k_mb = 0.02
alpha = 25
R = a
f_th = 0.0005

def create_adh(a, b, d, xc, yc, angle, max_num_adh, rng):

    num = rng.integers(0, max_num_adh, 1)[0]

    adh_arr0 = np.zeros((max_num_adh, 2))
    adh_arr = np.zeros((max_num_adh, 2))

    i = 0
    while (i < num):

        x = d + a*(2*rng.random(1)[0] - 1)
        y = b*(2*rng.random(1)[0] - 1)

        if incell(x, y, a, b, d):
            u = x*cos(angle) + y*sin(angle)
            v = -x*sin(angle) + y*cos(angle)

            adh_arr0[i] = np.array([u, v])
            adh_arr[i] = np.array([u, v])
            i = i + 1

    return adh_arr0, adh_arr

if __name__ == '__main__':

    a = 6
    b = 4
    d = 1
    xc = 0.0
    yc = 0.0
    angle = 0.0
    lamda = 0.4
    ks = 0.1
    max_num_adh = 64
    k_mf = 0.02
    k_mb = 0.02
    alpha = 25
    R = a
    f_th = 0.0005

    rng = np.random.default_rng(652)

    adh_arr0, adh_arr = create_adh(a, b, d, xc, yc, angle, max_num_adh, rng)

    # polygon
    phi = np.linspace(0, 2*pi, 1000)
    x_coords = ((a**2-d**2*np.sin(phi))**0.5 + d*np.cos(phi))*np.cos(phi)
    y_coords = b*np.sin(phi)

    coords = [[x, y] for x, y in list(zip(x_coords, y_coords))]
    polygon = shapely.geometry.Polygon(coords)

    fig, ax = plt.subplots(subplot_kw={'aspect': 'equal'})
    verts = np.array(polygon.exterior.coords.xy)
    patch = Polygon(verts.T, color = 'blue', alpha = 0.3)
    ax.add_patch(patch)

    ax.scatter(adh_arr[:, 0], adh_arr[:, 1], c = "red")

    ax.set_xlim(-2*a, xc+2*a)
    ax.set_ylim(-2*b, yc+2*b)

    plt.show()
