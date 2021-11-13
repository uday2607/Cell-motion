from matplotlib import pyplot as plt
from shapely.geometry.point import Point
from shapely import affinity
from scipy.spatial import Voronoi
from matplotlib.patches import Polygon
import numpy as np
import math

def voronoi_finite_polygons_2d(vor, radius=None):
    """
    Reconstruct infinite voronoi regions in a 2D diagram to finite
    regions.

    Parameters
    ----------
    vor : Voronoi
        Input diagram
    radius : float, optional
        Distance to 'points at infinity'.

    Returns
    -------
    regions : list of tuples
        Indices of vertices in each revised Voronoi regions.
    vertices : list of tuples
        Coordinates for revised Voronoi vertices. Same as coordinates
        of input vertices, with 'points at infinity' appended to the
        end.

    """

    if vor.points.shape[1] != 2:
        raise ValueError("Requires 2D input")

    new_regions = []
    new_vertices = vor.vertices.tolist()

    center = vor.points.mean(axis=0)
    if radius is None:
        radius = vor.points.ptp().max()

    # Construct a map containing all ridges for a given point
    all_ridges = {}
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))

    # Reconstruct infinite regions
    for p1, region in enumerate(vor.point_region):
        vertices = vor.regions[region]

        if all(v >= 0 for v in vertices):
            # finite region
            new_regions.append(vertices)
            continue

        # reconstruct a non-finite region
        ridges = all_ridges[p1]
        new_region = [v for v in vertices if v >= 0]

        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1
            if v1 >= 0:
                # finite ridge: already in the region
                continue

            # Compute the missing endpoint of an infinite ridge

            t = vor.points[p2] - vor.points[p1] # tangent
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])  # normal

            midpoint = vor.points[[p1, p2]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            far_point = vor.vertices[v2] + direction * radius

            new_region.append(len(new_vertices))
            new_vertices.append(far_point.tolist())

        # sort region counterclockwise
        vs = np.asarray([new_vertices[v] for v in new_region])
        c = vs.mean(axis=0)
        angles = np.arctan2(vs[:,1] - c[1], vs[:,0] - c[0])
        new_region = np.array(new_region)[np.argsort(angles)]

        # finish
        new_regions.append(new_region.tolist())

    return new_regions, np.asarray(new_vertices)

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

i = 0
while i < num:
    j = 1
    temp = np.multiply(np.array([X, Y, 360]), np.random.rand(3))
    temp_ell = create_ellipse((temp[0], temp[1]),(a,b),temp[2])

    # Current center
    theta = math.pi*temp[2]/180
    pot = temp[:2].copy() + ((4*a**2-b**2)**0.5)*np.array([math.cos(theta), math.sin(theta)])/2
    circ1 = Point(pot[:2]).buffer(a, resolution=1000)
    pot = temp[:2].copy() - ((4*a**2-b**2)**0.5)*np.array([math.cos(theta), math.sin(theta)])/2
    circ2 = Point(pot[:2]).buffer(a, resolution=1000)

    for cent in centers:

        theta = math.pi*cent[2]/180
        pot = cent[:2].copy() + ((4*a**2-b**2)**0.5)*np.array([math.cos(theta), math.sin(theta)])/2
        cir1 = Point(pot[:2]).buffer(a, resolution=1000)
        pot = cent[:2].copy() - ((4*a**2-b**2)**0.5)*np.array([math.cos(theta), math.sin(theta)])/2
        cir2 = Point(pot[:2]).buffer(a, resolution=1000)

        if ((circ1.intersects(cir1) or circ1.intersects(cir2)) or
            (circ2.intersects(cir1) or circ2.intersects(cir2))):
            j = 0
            break

    if (j == 1):
        centers.append(temp)
        ells.append(temp_ell)
        print(i)
        i = i+1

fig, ax = plt.subplots(subplot_kw={'aspect': 'equal'})
for e in ells:
    verts = np.array(e.exterior.coords.xy)
    patch = Polygon(verts.T, color = 'blue', alpha = 0.5)
    ax.add_patch(patch)

ax.set_xlim(-2*a, X+2*a)
ax.set_ylim(-2*b, Y+2*b)

## Cells in my model
points = []
for cent in centers:

    theta = math.pi*cent[2]/180
    temp_pts = []

    pot = cent[:2].copy() + ((4*a**2-b**2)**0.5)*np.array([math.cos(theta), math.sin(theta)])/2
    circ = Point(pot[:2]).buffer(a, resolution=1000)
    verts = np.array(circ.exterior.coords.xy)
    patch = Polygon(verts.T, color = 'red', alpha = 0.2)
    ax.add_patch(patch)
    points.append(pot)

    pot = cent[:2].copy() - ((4*a**2-b**2)**0.5)*np.array([math.cos(theta), math.sin(theta)])/2
    circ = Point(pot[:2]).buffer(a, resolution=1000)
    verts = np.array(circ.exterior.coords.xy)
    patch = Polygon(verts.T, color = 'red', alpha = 0.2)
    ax.add_patch(patch)
    points.append(pot)

ax.set_xlim(-2*a, X+2*a)
ax.set_ylim(-2*b, Y+2*b)


## Voronoi Teselltions
# compute Voronoi tesselation
points = np.asarray(points)
plt.plot(points[:,0],points[:,1],'w.')
plt.axis('square')
vor = Voronoi(points)
skmax = np.shape(vor.ridge_vertices)[0]
for sk in np.arange(skmax):
    a = vor.ridge_points[sk,0]
    b = vor.ridge_points[sk,1]
    c = min(a,b)
    d = abs(a-b)
    flag = 1
    if d == 1 and np.mod(c,2) == 0:
        flag = 0
    if flag == 1:
        s = vor.ridge_vertices[sk]
        s = np.asarray(s)
        if np.all(s>=0):
            plt.plot(vor.vertices[s,0],vor.vertices[s,1],'k-')
for k in np.arange(0,2*num,2):
    plt.plot([points[k,0],points[k+1,0]],[points[k,1],points[k+1,1]],'r-')
plt.show()
