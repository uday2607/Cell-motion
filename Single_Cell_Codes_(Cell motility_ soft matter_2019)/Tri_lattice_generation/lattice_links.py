import numpy as np
from matplotlib import pyplot as plt

import math, random
from hexalattice.hexalattice import create_hex_grid


# Create Hexagonal grid centers
L = 15
hex_centers_p, _ = create_hex_grid(nx=L+2, ny=L+2, do_plot=False)

# Create indices for the hex_centers
indices = np.arange((L+2)*(L+2))
indices_p = np.tile(np.arange((L+2)*(L+2)),((L+2)*(L+2), 1))

# Calculate the distance between each points
dist = np.around(np.linalg.norm(hex_centers_p - hex_centers_p[:,None], axis=-1), 2)

# Isolate the neraest neighbors
idx = np.where(dist == 1.0, 1, 0)
ind = indices[np.count_nonzero(idx == 1, axis=0) == 6]
idx[np.count_nonzero(idx == 1, axis=0) != 6] = np.zeros((L+2)*(L+2))

# Save the indices of the neighbors
nbs = indices_p[idx == 1].reshape((-1, 6))

# Plot
for i in range(len(ind)):
    plt.scatter(*hex_centers_p[ind[i]], c="red")
    for j in range(6):
        if random.random() != 0:
            xt, yt = zip(hex_centers_p[ind[i]], hex_centers_p[nbs[i]][j])
            plt.plot(xt, yt, c="blue")

plt.xlim([-L/2, L/2])
plt.ylim([-L/2, L/2])
plt.show()
