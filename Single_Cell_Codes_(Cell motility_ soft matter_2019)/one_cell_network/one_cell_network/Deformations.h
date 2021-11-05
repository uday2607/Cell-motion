void AffineShear(double *x, double ShearStrain) {
    //note for 2D case, new x coordinations also depends on y coordinations. (instead of z)
    int G1, N = Lx * Ly;
    for (int i = 0; i <= (Lx - 1); i++) {
        for (int j = 0; j <= (Ly - 1); j++) {
            G1 = 1 + j * Lx + i;

            // y = y0 * (1+eps)

            x[G1 + N + 0 * 2 * N] = (1 + eps) * x[G1 + N + 0 * 2 * N];   // y coordinates


            // x = x0 + y*strain

            x[G1 + 0 * 2 * N] = x[G1 + 0 * 2 * N] + x[G1 + N + 0 * 2 * N] * ShearStrain;    // x coordinates


        }
    }
    return;
}
