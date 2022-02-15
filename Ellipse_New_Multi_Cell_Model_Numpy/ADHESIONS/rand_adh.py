from Cell_funcs import *

"""New Adhesions"""
def random_adhesions(cells, cparams,
                    Nadh, k_plus, dt, rng):

    Adh = np.zeros((cells.shape[0]//3, Nadh, 2)) - 1e8
    Adh0 = np.zeros((cells.shape[0]//3, Nadh, 2)) - 1e8
    cAdh0 = np.zeros((cells.shape[0]//3, Nadh, 2)) - 1e8
    cp0 = np.zeros((cells.shape[0]//3, Nadh, 3)) - 1e8

    for ind in range(cells.shape[0]//3):
        rands = rng.uniform(0, 1, Nadh)

        #indices of adhesions to be formed
        adh_ind = np.arange(Nadh)[rands < k_plus*dt]

        for ad in adh_ind:

            # Find a point from the grid
            found = 0

            while found == 0:

                phi = rng.uniform(0, 2*np.pi)
                rho = rng.uniform(0.0, 1.0)
                x_ = cparams[4*ind]*np.sqrt(rho)*np.cos(phi)
                y_ = cparams[4*ind+1]*np.sqrt(rho)*np.sin(phi)
                x = x_*cos(cparams[4*ind+2]) - y_*sin(cparams[4*ind+2])
                y = x_*sin(cparams[4*ind+2]) + y_*cos(cparams[4*ind+2])
                x = x_+cells[3*ind]
                y = y_+cells[3*ind+1]

                if (in_ellipse(x, y, cells[3*ind], cells[3*ind+1],
                    cparams[4*ind], cparams[4*ind+1], cparams[4*ind+2])):

                    # Check if there is an adhesion point from this site
                    if not np.any(np.logical_and(Adh[ind, :, 0] == x,
                            Adh[ind, :, 1] == y)):
                        found = 1

            #Store adhesions
            #From ground reference
            Adh0[ind, ad, 0] = x
            Adh0[ind, ad, 1] = y
            Adh[ind, ad, 0] = x - cells[3*ind]
            Adh[ind, ad, 1] = y - cells[3*ind+1]
            cAdh0[ind, ad] = Adh[ind, ad].copy()
            cp0[ind, ad, 0] = cparams[4*ind]
            cp0[ind, ad, 1] = cparams[4*ind+2]
            cp0[ind, ad, 2] = 0

    return Adh0, Adh, cAdh0, cp0

def one_cell_random_adh(cells, ind, cparams, Adh, Adh0, cAdh0, cp0,
                    Nadh, k_plus, dt, flag, rng):

    # Add adhesions to existing ones
    if flag == 0:
        #Indices of adhesion sites not bonded yet
        adh_ind = np.arange(Adh.shape[1])[np.logical_and(
            Adh[ind, :, 0] == -1e8, Adh[ind, :, 1] == -1e8)]

        #indices of adhesions to be formed
        rands = rng.uniform(0, 1, adh_ind.size)
        adh_ind = adh_ind[rands < k_plus*dt]

        for ad in adh_ind:

            # Find a point from the grid
            found = 0

            while found == 0:
                if flag == 0:
                    phi = rng.uniform(0, 2*np.pi)
                    rho = rng.uniform(0.0, 1.0)
                else:
                    phi = rng.uniform(0, np.pi) - np.pi/2
                    rho = rng.uniform(0.0, 1.0)

                x_ = cparams[4*ind]*sqrt(rho)*cos(phi)
                y_ = cparams[4*ind+1]*sqrt(rho)*sin(phi)
                x = x_*cos(cparams[4*ind+2]) - y_*sin(cparams[4*ind+2])
                y = x_*sin(cparams[4*ind+2]) + y_*cos(cparams[4*ind+2])
                x = x_+cells[3*ind]
                y = y_+cells[3*ind+1]

                if (in_ellipse(x, y, cells[3*ind], cells[3*ind+1],
                    cparams[4*ind], cparams[4*ind+1], cparams[4*ind+2])):

                    # Check if there is an adhesion point from this site
                    if not np.any(np.logical_and(Adh[ind, :, 0] == x,
                            Adh[ind, :, 1] == y)):
                        found = 1

            #Store adhesions
            #From ground reference
            Adh0[ind, ad, 0] = x
            Adh0[ind, ad, 1] = y
            Adh[ind, ad, 0] = x - cells[3*ind]
            Adh[ind, ad, 1] = y - cells[3*ind+1]
            cAdh0[ind, ad] = Adh[ind, ad].copy()
            cp0[ind, ad, 0] = cparams[4*ind]
            cp0[ind, ad, 1] = cparams[4*ind+2]
            cp0[ind, ad, 2] = 0

    # completely new adhesions
    elif flag == 1:

        Adh[ind] = np.zeros((Nadh, 2)) - 1e8
        Adh0[ind] = np.zeros((Nadh, 2)) - 1e8
        cAdh0[ind] = np.zeros((Nadh, 2)) - 1e8
        cp0[ind] = np.zeros((Nadh, 3)) - 1e8

        #Indices of adhesion sites not bonded yet
        adh_ind = np.arange(Nadh)

        #indices of adhesions to be formed
        rands = rng.uniform(0, 1, Nadh)
        adh_ind = adh_ind[rands < k_plus*dt]

        for ad in adh_ind:

            # Find a point from the grid
            found = 0

            while found == 0:
                if flag == 0:
                    phi = rng.uniform(0, 2*np.pi)
                    rho = rng.uniform(0.0, 1.0)
                else:
                    phi = rng.uniform(0, np.pi) - np.pi/2
                    rho = rng.uniform(0.0, 1.0)

                x_ = cparams[4*ind]*sqrt(rho)*cos(phi)
                y_ = cparams[4*ind+1]*sqrt(rho)*sin(phi)
                x = x_*cos(cparams[4*ind+2]) - y_*sin(cparams[4*ind+2])
                y = x_*sin(cparams[4*ind+2]) + y_*cos(cparams[4*ind+2])
                x = x_+cells[3*ind]
                y = y_+cells[3*ind+1]

                if (in_ellipse(x, y, cells[3*ind], cells[3*ind+1],
                    cparams[4*ind], cparams[4*ind+1], cparams[4*ind+2])):

                    # Check if there is an adhesion point from this site
                    if not np.any(np.logical_and(Adh[ind, :, 0] == x,
                            Adh[ind, :, 1] == y)):
                        found = 1

            #Store adhesions
            #From ground reference
            Adh0[ind, ad, 0] = x
            Adh0[ind, ad, 1] = y
            Adh[ind, ad, 0] = x - cells[3*ind]
            Adh[ind, ad, 1] = y - cells[3*ind+1]
            cAdh0[ind, ad] = Adh[ind, ad].copy()
            cp0[ind, ad, 0] = cparams[4*ind]
            cp0[ind, ad, 1] = cparams[4*ind+2]
            cp0[ind, ad, 2] = 0
    
    return Adh0[ind], Adh[ind], cAdh0[ind], cp0[ind]
""""""