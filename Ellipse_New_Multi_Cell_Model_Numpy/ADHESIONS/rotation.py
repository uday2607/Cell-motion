from Cell_funcs import *

"""Rotation and shift of Adhesion sites"""
@nb.jit(nopython = True, nogil = True)
def rotation_and_shift(cells, cell_p, cparams, Adh):

    for i in range(cells.shape[0]//3):
        ind = np.arange(Adh.shape[1])[np.logical_and(
            Adh[i, :, 0] != -1e8, Adh[i, :, 1] != -1e8)]

        if ind.shape[0] != 0:

            #rotate all the adhesions
            Adh_ = Adh[i].copy()
            x, y = Adh_[ind, 0], Adh_[ind, 1]
            dtheta = cells[3*i+2]
            Adh_[ind, 0] = np.cos(dtheta)*x - np.sin(dtheta)*y
            Adh_[ind, 1] = np.sin(dtheta)*x + np.cos(dtheta)*y
            Adh[i] = Adh_

        #rotate theta
        cparams[4*i+2] += cells[3*i+2]
        cells[3*i+2] = 0

    return cells, cparams, Adh
""""""