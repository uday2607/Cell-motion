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

## Normal Ellipses
X = 100
Y = 100
num = 50
a = 4
b = 3

coord = np.multiply(np.array([X, Y, 360]), np.random.rand(3))
centers = [coord]
ells = [create_ellipse(coord[:2],(a,b),coord[2])]
line_segments = []

i = 0
while i < num:
    j = 1
    temp = np.multiply(np.array([X, Y, 360]), np.random.rand(3))
    temp_ell = create_ellipse((temp[0], temp[1]),(a,b),temp[2])

    for e in ells:
        if temp_ell.boundary.intersects(e.boundary):
            points = list(temp_ell.boundary.intersection(e.boundary).geoms)
            lsg = [LineString([p1, p2]) for p1, p2 in list(
                        combinations(points, 2))]
            line_segments = line_segments+lsg

    centers.append(temp)
    ells.append(temp_ell)
    print(i)
    i = i+1


fig, ax = plt.subplots(subplot_kw={'aspect': 'equal'})
ells = unary_union(ells)
for e in polygonize(ells):
    verts = np.array(e.exterior.coords.xy)
    patch = Polygon(verts.T, color = 'blue', alpha = 0.5)
    ax.add_patch(patch)

for lg in line_segments:
    ax.plot(*lg.xy, color = "red")

ax.set_xlim(-2*a, X+2*a)
ax.set_ylim(-2*b, Y+2*b)

plt.show()
