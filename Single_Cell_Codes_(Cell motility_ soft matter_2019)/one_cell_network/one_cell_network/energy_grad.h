double AdEnergy(double *x) {
    int num = 1;
    int N = Lx * Ly;

    double xc_new = x[2 * N + 1];  //nvariables = 1 + 2*Lx*Ly+3;
    double yc_new = x[2 * N + 2];
    double theta = x[2 * N + 3];

    x[1] = 0;     //fix one corner of the network to prevent possible shift
    x[1 + N] = 0;

    double E = 0;

    while (num <= totalnum) {
        if (ifremove[num] == 1) {
            num++;
            continue;
        }   //the adhesion site has been removed

        double xx = Adhesion[num][0];
        double yy = Adhesion[num][1];

        int ix, iy;
        ix = AdhesionIndex[num][0];
        iy = AdhesionIndex[num][1];

        double x0 = x[ix] + Adhesion_shift[num][0];    //now x0, y0 is no longer fixed
        double y0 = x[iy] + Adhesion_shift[num][1];    //shift x0, y0 due to periodic boundary condition

        double x_new = cos(theta) * xx - sin(theta) * yy;
        double y_new = sin(theta) * xx + cos(theta) * yy;

        double fx = -(x_new + xc_new - x0);
        double fy = -(y_new + yc_new - y0);

//        if (fabs(fx)>0.1||fabs(fy)>0.1){
//            cout<<num<<endl;
//            cout<<fx<<" "<<fy<<endl;
//            cout<<"----------"<<endl;
//        }

        E = E + 0.5 * ks * fx * fx + 0.5 * ks * fy * fy;

        num++;

    }

    return E;

}

/*
double AdEnergy2(double *x) {
    int num = 1;
    int N = Lx * Ly;

    double xc_new = x[2 * N + 4];  //nvariables = 1 + 2*Lx*Ly+3;
    double yc_new = x[2 * N + 5];
    double theta = x[2 * N + 6];

    x[1] = 0;     //fix one corner of the network to prevent possible shift
    x[1 + N] = 0;

    double E = 0;

    while (num <= totalnum2) {
        if (ifremove2[num] == 1) {
            num++;
            continue;
        }   //the adhesion site has been removed

        double xx = Adhesion2[num][0];
        double yy = Adhesion2[num][1];

        int ix, iy;
        ix = AdhesionIndex2[num][0];
        iy = AdhesionIndex2[num][1];

        double x0 = x[ix] + Adhesion_shift2[num][0];    //now x0, y0 is no longer fixed
        double y0 = x[iy] + Adhesion_shift2[num][1];    //shift x0, y0 due to periodic boundary condition

        double x_new = cos(theta) * xx - sin(theta) * yy;
        double y_new = sin(theta) * xx + cos(theta) * yy;

        double fx = -(x_new + xc_new - x0);
        double fy = -(y_new + yc_new - y0);

//        cout<<"fx="<<fx<<endl;
//        cout<<"fy="<<fy<<endl;

//        if (fabs(fx)>0.1||fabs(fy)>0.1){
//            cout<<num<<endl;
//            cout<<fx<<" "<<fy<<endl;
//            cout<<"----------"<<endl;
//        }

        E = E + 0.5 * ks * fx * fx + 0.5 * ks * fy * fy;

        num++;

    }

    return E;

}
*/

void Adgrad(double *x, double *g) {

    int num = 1;
    int N = Lx * Ly;

    double xc_new = x[2 * N + 1];  //nvariables = 1 + 2*Lx*Ly+3;
    double yc_new = x[2 * N + 2];
    double theta = x[2 * N + 3];

    x[1] = 0;     //fix one corner of the network to prevent possible shift
    x[1 + N] = 0;

    while (num <= totalnum) {

        if (ifremove[num] == 1) {
            num++;
            continue;
        }

        double xx = Adhesion[num][0];
        double yy = Adhesion[num][1];

        double x0 = x[AdhesionIndex[num][0]] + Adhesion_shift[num][0];    //now x0, y0 is no longer fixed
        double y0 = x[AdhesionIndex[num][1]] + Adhesion_shift[num][1];

        double x_new = cos(theta) * xx - sin(theta) * yy;
        double y_new = sin(theta) * xx + cos(theta) * yy;

        double fx = -(x_new + xc_new - x0);
        double fy = -(y_new + yc_new - y0);

        double torque = -(xc_new - x0) * y_new + (yc_new - y0) * x_new;

//        E=E+0.5*fx*fx+0.5*fy*fy;

        g[AdhesionIndex[num][0]] = g[AdhesionIndex[num][0]] + ks * fx;
        g[AdhesionIndex[num][1]] = g[AdhesionIndex[num][1]] + ks * fy;

//        if (fabs(fx)>0.1||fabs(fy)>0.1||fabs(torque)>0.1){
//            cout<<num<<endl;
//            cout<<ks*fx<<" "<<ks*fy<<" "<<ks*torque<<endl;
//            cout<<"----------"<<endl;
//        }

        g[2 * N + 1] = g[2 * N + 1] - ks * fx;
        g[2 * N + 2] = g[2 * N + 2] - ks * fy;
        g[2 * N + 3] = g[2 * N + 3] + ks * torque;

        num++;
    }

    g[1] = 0.0;
    g[1 + N] = 0.0;
    return;
}


/*
void Adgrad2(double *x, double *g) {

    int num = 1;
    int N = Lx * Ly;

    double xc_new = x[2 * N + 4];  //nvariables = 1 + 2*Lx*Ly+3;
    double yc_new = x[2 * N + 5];
    double theta = x[2 * N + 6];

    x[1] = 0;     //fix one corner of the network to prevent possible shift
    x[1 + N] = 0;

    while (num <= totalnum2) {

        if (ifremove2[num] == 1) {
            num++;
            continue;
        }

        double xx = Adhesion2[num][0];
        double yy = Adhesion2[num][1];

        double x0 = x[AdhesionIndex2[num][0]] + Adhesion_shift2[num][0];    //now x0, y0 is no longer fixed
        double y0 = x[AdhesionIndex2[num][1]] + Adhesion_shift2[num][1];

        double x_new = cos(theta) * xx - sin(theta) * yy;
        double y_new = sin(theta) * xx + cos(theta) * yy;

        double fx = -(x_new + xc_new - x0);
        double fy = -(y_new + yc_new - y0);

        double torque = -(xc_new - x0) * y_new + (yc_new - y0) * x_new;

//        E=E+0.5*fx*fx+0.5*fy*fy;

        g[AdhesionIndex2[num][0]] = g[AdhesionIndex2[num][0]] + ks * fx;
        g[AdhesionIndex2[num][1]] = g[AdhesionIndex2[num][1]] + ks * fy;

//        if (fabs(fx)>0.1||fabs(fy)>0.1||fabs(torque)>0.1){
//            cout<<num<<endl;
//            cout<<ks*fx<<" "<<ks*fy<<" "<<ks*torque<<endl;
//            cout<<"----------"<<endl;
//        }

        g[2 * N + 4] = g[2 * N + 4] - ks * fx;
        g[2 * N + 5] = g[2 * N + 5] - ks * fy;
        g[2 * N + 6] = g[2 * N + 6] + ks * torque;

        num++;
    }

    g[1] = 0.0;
    g[1 + N] = 0.0;
    return;
}
*/


#define tiny 1e-16

double powMy(double A, int n) {
    double temp = 1;
    for (int i = 1; i <= n; i++) {
        temp = temp * A;
    }
    return temp;
}

double EnergyStretch(double *x) {
    int i, j, con, G1, G2, N = Lx * Ly;
    double Energy = 0, dx, dy, dz, dis;

    x[1] = 0;     //fix one corner of the network to prevent possible shift
    x[1 + N] = 0;

    for (i = 0; i <= (Lx - 1); i++) {
        for (j = 0; j <= (Ly - 1); j++) {

            G1 = 1 + j * Lx + i;
            for (con = 0; con <= 2; con++) {
                G2 = Gcon[con][G1];
                if (G2 != -1) {
                    //-----beautiful code for phantomalization---------------

                    dx = x[G2] - x[G1];
                    dy = x[G2 + N] - x[G1 + N];
                    //-----beautiful code for phantomalization---------------

                    //for con==1,2,3 all the following codes are valid, since they only change one index at a time

                    if (BoundaryX[con][G1]) {
                        dx = dx + Lx * (con == 0 ? -1.0 : 1.0);

                    }
                    if (BoundaryY[con][G1]) {
                        dy = dy + h;   //h=Ly*0.86602540378443864676*(1+eps);

                        dx = dx + 0.5 * Ly;
//                        dx=dx+Ly*0.866025403784438548717*ShearStrain;   //lee-edwards boundary condition
                    }

                    dis = sqrt(dx * dx + dy * dy);
                    Energy += Kcon[con][G1] * powMy(dis - 1.0, 2) / 2.0;
                }
            }

        }
    }

    return Energy;
}


void StretchEnergyGrad(double *x, double *Grad) {

    double dx, dy, temp;
    int i, j, con, G1, G2, G3, N = Lx * Ly;
    double x1, x2, x3, y1, y2, y3, dis;

    x[1] = 0;     //fix one corner of the network to prevent possible shift
    x[1 + N] = 0;

    for (i = 0; i < nvariables; i++)   //the last three elements are xc, yc and theta
    {
        Grad[i] = 0.0;
    }

    for (i = 0; i <= (Lx - 1); i++) {
        for (j = 0; j <= (Ly - 1); j++) {
            G1 = 1 + j * Lx + i;
            for (con = 0; con <= 2; con++) {
                G2 = Gcon[con][G1];
                if (G2 != -1) {
                    dx = x[G2] - x[G1];
                    dy = x[G2 + N] - x[G1 + N];

                    if (BoundaryX[con][G1]) {
                        dx = dx + Lx * (con == 0 ? -1.0 : 1.0);
                    }
                    if (BoundaryY[con][G1]) {
                        dy = dy + h; //h=Ly*0.86602540378443864676*(1+eps);

                        dx = dx + 0.5 * Ly;

//                        dx=dx+Ly*0.866025403784438548717*ShearStrain;   //lee-edwards boundary condition
                    }

                    dis = sqrt(dx * dx + dy * dy);

                    temp = Kcon[con][G1] * (dis - 1.0) / (dis + 1e-16);

                    Grad[G1 + 0 * N] += -dx * temp;
                    Grad[G1 + 1 * N] += -dy * temp;

                    Grad[G2 + 0 * N] += dx * temp;
                    Grad[G2 + 1 * N] += dy * temp;

                    G3 = Gcon[con][G2];
                    if (G3 != -1) {
                        x1 = x[G1];
                        y1 = x[G1 + N];
                        x2 = x[G2];
                        y2 = x[G2 + N];
                        x3 = x[G3];
                        y3 = x[G3 + N];

                        if (BoundaryX[con][G1])   //G1:x2, x3  G2: x3
                        {
                            x2 += Lx * (con == 0 ? -1 : 1);
                            x3 += Lx * (con == 0 ? -1 : 1);
                        } else if (BoundaryX[con][G2]) {
                            x3 += Lx * (con == 0 ? -1 : 1);
                        }

                        if (BoundaryY[con][G1]) {
                            y2 += h;  //h=Ly*0.86602540378443864676*(1+eps);
                            y3 += h;

                            x2 += 0.5 * Ly;
                            x3 += 0.5 * Ly;
                        } else if (BoundaryY[con][G2]) {
                            y3 += h;  //h=Ly*0.86602540378443864676*(1+eps);
                            x3 += 0.5 * Ly;
                        }

                        double xji, yji, xjk, yjk;

                        xji = x2 - x1;
                        yji = y2 - y1;

                        xjk = x2 - x3;
                        yjk = y2 - y3;

                        double lji, ljk, dot;

                        dot = xji * xjk + yji * yjk;

                        lji = sqrt(xji * xji + yji * yji) + tiny;
                        ljk = sqrt(xjk * xjk + yjk * yjk) + tiny;

                        double temp = dot / (lji * ljk);

                        if (temp < -1.0) {
                            temp = -1.0;
                        } else if (temp > 1.0) {
                            temp = 1.0;
                        }

                        double angle = acos(temp) - M_PI, common_factor;

                        common_factor = kappa * angle * (-1 / (sqrt(1 - temp * temp) + tiny));

                        double c1, c2, c3, c4; //I need to check the gradient here.

                        c1 = xjk / (lji * ljk) - (xji * dot) / (ljk * lji * lji * lji);
                        c2 = xji / (lji * ljk) - (xjk * dot) / (lji * ljk * ljk * ljk);

                        c3 = yjk / (lji * ljk) - (yji * dot) / (ljk * lji * lji * lji);  //change x->y in c1;
                        c4 = yji / (lji * ljk) - (yjk * dot) / (lji * ljk * ljk * ljk);

                        Grad[G1 + 0 * N] += common_factor * (-c1);     //g[x1]
                        Grad[G2 + 0 * N] += common_factor * (c1 + c2);  //g[x2]
                        Grad[G3 + 0 * N] += common_factor * (-c2);     //g[x3]


                        Grad[G1 + 1 * N] += common_factor * (-c3);     //g[y1]
                        Grad[G2 + 1 * N] += common_factor * (c3 + c4);  //g[y2]
                        Grad[G3 + 1 * N] += common_factor * (-c4);     //g[y3]

                    }

                }
            }
        }
    }

    Grad[1] = 0.0;
    Grad[1 + N] = 0.0;

    return;
}

double cal_phi(double x_i, double y_i, double x_j, double y_j, double x_k, double y_k) {

    double temp, y;

    temp = ((x_j - x_i) * (x_j - x_k) + (y_j - y_i) * (y_j - y_k)) /
           ((sqrt((x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i)) + tiny) *
            (sqrt((x_j - x_k) * (x_j - x_k) + (y_j - y_k) * (y_j - y_k)) + tiny));

    if (temp < -1.0) {
        temp = -1.0;
    } else if (temp > 1.0) {
        temp = 1.0;
    }

    y = acos(temp) - M_PI;

    return y;

}

double EnergyBend(double x[]) {
    double Energy = 0, x1, x2, x3, y1, y2, y3;
    int i, j, con, G1, G2, G3, N = Lx * Ly;

    for (i = 0; i <= (Lx - 1); i++) {
        for (j = 0; j <= (Ly - 1); j++) {
            G1 = 1 + j * Lx + i;
            for (con = 0; con <= 2; con++) {
                G2 = Gcon[con][G1];
                if (G2 != -1) {
                    G3 = Gcon[con][G2];
                    if (G3 != -1) {
                        x1 = x[G1];
                        y1 = x[G1 + N];
                        x2 = x[G2];
                        y2 = x[G2 + N];
                        x3 = x[G3];
                        y3 = x[G3 + N];

                        if (BoundaryX[con][G1])   //G1:x2, x3  G2: x3
                        {
                            x2 += Lx * (con == 0 ? -1 : 1);
                            x3 += Lx * (con == 0 ? -1 : 1);
                        } else if (BoundaryX[con][G2]) {
                            x3 += Lx * (con == 0 ? -1 : 1);
                        }

                        if (BoundaryY[con][G1]) {
                            y2 += h;  //h=Ly*0.86602540378443864676*(1+eps);
                            y3 += h;

                            x2 += 0.5 * Ly;
                            x3 += 0.5 * Ly;
                        } else if (BoundaryY[con][G2]) {
                            y3 += h;
                            x3 += 0.5 * Ly;
                        }

                        double angle = cal_phi(x1, y1, x2, y2, x3, y3);

                        Energy = Energy + 0.5 * angle * angle;
                    }
                }
            }

        }
    }
    return Energy;
}

double energy(double *x) {
    double E;

    E = AdEnergy(x) /*+ AdEnergy2(x)*/ + EnergyStretch(x) + kappa * EnergyBend(x);

    return E;
}

void grad(double *x, double *g) {

    StretchEnergyGrad(x, g);  //contains initialization of g
    Adgrad(x, g);
//    Adgrad2(x, g);
}

void Numerical_EnergyGrad(double *x, double *Grad) {
    double x0[4];
    double energy1, energy2;

    for (int i = 0; i < 4; i++) {
        x0[i] = x[i];
    }

    for (int j = 0; j < 4; j++) {

        x0[j] = x[j] + 1e-4;

        energy1 = energy(x0);

        x0[j] = x[j] - 1e-4;

        energy2 = energy(x0);

        Grad[j] = (energy1 - energy2) / 2e-4;

        x0[j] = x[j];
    }

}