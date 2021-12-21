import cell
import adhesions
import numpy as np

if __name__ == '__main__':
    a = 3
    b = 2
    L = 40
    Num = 4
    Nadh = 32
    k_plus = 0.4
    dt = 1
    lamda = 0.4
    tau = 30
    T_S = 10
    rng = np.random.default_rng()

    cells, cparams, Ovlaps = cell.random_cells(L, a, b, Num, rng)
    Adh0, Adh = adhesions.random_adhesions(L, a, b, cells, Nadh,
                                           k_plus, dt, rng)
    adhesions.contraction(cells, Ovlaps, Adh, lamda, tau, dt, T_S, a, b)

    print(Adh0)
