from matplotlib import pyplot as plt
import numpy as np
import math


a = 4
b = 2
L = 1000
Num = 10
Nadh = 64
k_plus = 0.025
lamda = 0.2
tau = 30
dt = 60/tau
T_S_c = 20000
T_S_p = 20000
k_s = 50.0
k_m = 1
E_const = np.zeros((3, 3))
E_const[2, 2] = -0.025
E_const[1, 1] = 2.0
E_const[1, 0] = 5000.0
E_const[0, 1] = 5000.0
E_const[0, 0] = 10000.0
fThreshold = 40.0/k_s
k_b = 0.025
k_f = 0.0025
alpha = 25
a_min = a
for i in range(tau):
    a_min -= a_min*lamda*dt/tau
rng = np.random.default_rng()
TIME = 600

with open("data.npy", "rb") as f:
    CELLS = np.load(f)
    CPARAMS = np.load(f)
    ADH = np.load(f)
    ADH0 = np.load(f)
    CADH0 = np.load(f)
    CP0 = np.load(f)


for t in range(TIME):

    X_v = []
    Y_v = []
    U_v = []
    V_v = []

    cells_x, cells_y = CELLS[t, ::3], CELLS[t, 1::3]

    x_min = int(np.min(cells_x))-5*a
    x_max = int(np.max(cells_x))+5*a
    y_min = int(np.min(cells_y))-5*a
    y_max = int(np.max(cells_y))+5*a

    steps = 5
    X, Y = np.meshgrid(np.arange(x_min, x_max, step=1/steps),
            np.arange(y_min, y_max, step=1/steps))
    Z = np.zeros((X.shape[1], Y.shape[0]))

    fig, ax = plt.subplots()

    Adh = ADH[t]
    Adh0 = ADH0[t]
    for num in range(cells_x.shape[0]):
        ind = np.arange(Adh.shape[1])[np.logical_and(
                Adh[num, :, 1] != -1e8, Adh[num, :, 0] != -1e8)]

        for i in ind:
            x_ad, y_ad = (steps*Adh0[num, i, 0])//steps, (steps*Adh0[num, i, 1])//steps
            print(x_ad, Adh[num, i, 0] + cells_x[num])
            X_v.append(x_ad)
            Y_v.append(y_ad)
            U_v.append(Adh[num, i, 0] + cells_x[num])
            V_v.append(Adh[num, i, 1] + cells_y[num])

            f = k_s*(np.sqrt((Adh[num, i, 0] + cells_x[num] - x_ad)**2
                    + np.sqrt((Adh[num, i, 1] + cells_y[num] - y_ad)**2)))

            z_ind_x, z_ind_y = steps*int(x_ad - x_min), steps*int(y_ad - y_min)
            Z[z_ind_x, z_ind_y] += f

    z_min, z_max = -np.abs(Z).max(), np.abs(Z).max()
    #c = ax.pcolormesh(X, Y, Z.T, vmin=z_min, vmax=z_max)
    ax.quiver(*np.array([X_v, Y_v]), np.array(U_v), np.array(V_v), units='xy',
                angles="xy", scale=1, linewidth=.5)
    #ax.set_title('pcolormesh')
    # set the limits of the plot to the limits of the data
    ax.axis([X.min(), X.max(), Y.min(), Y.max()])
    #fig.colorbar(c, ax=ax)
    ax.set_aspect(1./ax.get_data_ratio())
    plt.show()
