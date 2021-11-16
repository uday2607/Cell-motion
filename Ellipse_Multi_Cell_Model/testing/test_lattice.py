from matplotlib import pyplot as plt
from shapely.geometry.point import Point
from shapely import affinity
from scipy.spatial import Voronoi
from matplotlib.patches import Polygon
import numpy as np
import math

# Force the cell into the network
def inNet(L, a, b, x, y, theta):

    flag = 1

    # Compute extremes of the ellipse
    p1 = (x - a*math.cos(theta*math.pi/180),
          y - a*math.sin(theta*math.pi/180))
    p2 = (x + a*math.cos(theta*math.pi/180),
          y + a*math.sin(theta*math.pi/180))
    p3 = (x - b*math.sin(theta*math.pi/180),
          y + b*math.cos(theta*math.pi/180))
    p4 = (x + b*math.sin(theta*math.pi/180),
          y - b*math.cos(theta*math.pi/180))

    # Check if the cell is in the lattice
    # by checking if the cells are on the same side
    # of the lattice line segments
    if not (((p1[1] < math.sqrt(3)*p1[0]) and
        (p2[1] < math.sqrt(3)*p2[0]) and
        (p3[1] < math.sqrt(3)*p3[0]) and
        (p4[1] < math.sqrt(3)*p4[0]))):

        # push x-coordinate by a
        x = x + a
        flag = 0

    if not (((p1[1] > math.sqrt(3)*p1[0] - 2*L) and
        (p2[1] > math.sqrt(3)*p2[0] - 2*L) and
        (p3[1] > math.sqrt(3)*p3[0] - 2*L) and
        (p4[1] > math.sqrt(3)*p4[0] - 2*L))):

        # push x-coordinate by a
        x = x - a
        flag = 0

    if not ((p1[1] > 0) and (p2[1] > 0) and
        (p3[1] > 0) and (p4[1] > 0)):

        #push y-coordinate by b
        y = y + b
        flag = 0

    if not ((p1[1] < L) and (p2[1] < L) and
        (p3[1] < L) and (p4[1] < L)):

        #push y-coordinate by b
        y = y - b
        flag = 0

    return x, y, flag


def inellipse(x, y, theta, a, b, xi, yi):

    return (((x-xi)*math.cos(theta*math.pi/180) + (y-yi)*math.sin(theta*math.pi/180))**2/a**2 +
            ((x-xi)*math.sin(theta*math.pi/180) - (y-yi)*math.cos(theta*math.pi/180))**2/b**2 <= 1)

def noCollision(cell_coords, a, b, x, y, theta):

    NUM = 5000
    phi = np.linspace(0, 2*math.pi, NUM)
    ell1 = np.zeros((NUM, 2))
    ell1[:, 0] = x + a*math.cos(theta*math.pi/180)*np.cos(phi) - b*math.sin(theta*math.pi/180)*np.sin(phi)
    ell1[:, 1] = y + a*math.cos(theta*math.pi/180)*np.sin(phi) + b*math.sin(theta*math.pi/180)*np.cos(phi)

    for coord in cell_coords:
        x1, y1, theta1 = coord
        ell2 = np.zeros((NUM, 2))
        ell2[:, 0] = x1 + a*math.cos(theta1*math.pi/180)*np.cos(phi) - b*math.sin(theta1*math.pi/180)*np.sin(phi)
        ell2[:, 1] = y1 + a*math.cos(theta1*math.pi/180)*np.sin(phi) + b*math.sin(theta1*math.pi/180)*np.cos(phi)

        for p in ell1:
            xi, yi = p
            if inellipse(x1, y1, theta1, a, b, xi, yi):
                return False

        for p in ell2:
            xi, yi = p
            if inellipse(x, y, theta, a, b, xi, yi):
                return False

    return True

def create_ellipse(center, lengths, angle=0):
    """
    create a shapely ellipse. adapted from
    https://gis.stackexchange.com/a/243462
    """
    circ = Point(center).buffer(1, resolution=1000)
    ell = affinity.scale(circ, int(lengths[0]), int(lengths[1]))
    ellr = affinity.rotate(ell, angle)
    return ellr

L = 30
a = 3
b = 2
N = 10

fig, ax = plt.subplots(subplot_kw={'aspect': 'equal'})

points = np.zeros((L, L, 2))
for i in range(L):
    for j in range(L):
        points[i][j][0] = (2*i + j)*np.sqrt(1/3)
        points[i][j][1] = j

        ax.scatter(points[i][j][0], points[i][j][1], c="blue")

cell_coords = []
i = 0
while i < N:
    print(i)
    x, y, theta = np.array([np.random.uniform(a, L-a),
                  np.random.uniform(b, L-b), np.random.uniform(0, 360)])
    x = (2*x + y)/math.sqrt(3)

    # push the cell into the net
    flag = 0
    while (flag == 0):
        x, y, flag = inNet(L, a, b, x, y, theta)
    if noCollision(cell_coords, a, b, x, y, theta):
        cell_coords.append((x, y, theta))
        i += 1


ells = [create_ellipse(coord[:2],(a,b),coord[2]) for coord in cell_coords]

for e in ells:
    verts = np.array(e.exterior.coords.xy)
    patch = Polygon(verts.T, color = 'red', alpha = 0.5)
    ax.add_patch(patch)

plt.show()
