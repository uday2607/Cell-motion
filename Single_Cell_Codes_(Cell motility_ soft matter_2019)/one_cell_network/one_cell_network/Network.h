void Boundary(int &i2, int &j2)    //gives the corresponding points according to periodic boundary condition
{
    if (j2 < 0)
        j2 = j2 + Ly;
    else if (j2 >= Ly)
        j2 = j2 - Ly;

    if (i2 < 0)
        i2 = i2 + Lx;
    else if (i2 >= Lx)
        i2 = i2 - Lx;

    return;
}


void ConnectionsBuildFull(int **Gcon) {
    int i, j, i2, j2, G1, it, jt;

    // Bulk connection
    for (i = 0; i <= (Lx - 1); i++) {
        for (j = 0; j <= (Ly - 1); j++) {
            G1 = 1 + j * Lx + i;    //the index of each node is 1+Lx*Ly*z+Lx*y+x;
            //IN THE XY PLAIN
            // 0

            //--------- o is current node, v is (i2, j2, k2)----------
            //    v * * *
            //   * o * *
            //  * * * *
            //--------------------------------------------------------

            i2 = i - 1;
            j2 = j + 1;
            it = i2;
            jt = j2;
            Boundary(i2, j2);
            if (i2 != it) { BoundaryX[0][G1] = true; }
            else { BoundaryX[0][G1] = false; }
            if (j2 != jt) { BoundaryY[0][G1] = true; }
            else { BoundaryY[0][G1] = false; }

            Gcon[0][G1] = 1 + j2 * Lx + i2;    // Gcon[i][0] are meaningless
//            if (fRand(0,1)<pbond) Gcon[0][G1]=1+j2*Lx+i2;
//            else Gcon[0][G1]=-1;
//           if (j<=Ly/2-5) Kcon[0][G1]=1e-1;
//           else if (j>Ly/2+5) Kcon[0][G1]=1e1;
//           else{
//               Kcon[0][G1]=1e-1*pow(10, (j-Ly/2+5)*0.2);
//            }
            Kcon[0][G1] = 1;
            // 1

            //--------- o is current node, v is (i2, j2, k2)----------
            //    * v * *
            //   * o * *
            //  * * * *
            //--------------------------------------------------------

            i2 = i;
            j2 = j + 1;
            it = i2;
            jt = j2;
            Boundary(i2, j2);
            if (i2 != it) { BoundaryX[1][G1] = true; }
            else { BoundaryX[1][G1] = false; }
            if (j2 != jt) { BoundaryY[1][G1] = true; }
            else { BoundaryY[1][G1] = false; }
            Gcon[1][G1] = 1 + j2 * Lx + i2;
//            if (fRand(0,1)<pbond) Gcon[1][G1]=1+j2*Lx+i2;
//            else Gcon[1][G1]=-1;
//           if (j<=Ly/2-5) Kcon[1][G1]=1e-1;
//           else if (j>Ly/2+5) Kcon[1][G1]=1e1;
//           else{
//               Kcon[1][G1]=1e-1*pow(10, (j-Ly/2+5)*0.2);
//           }
            Kcon[1][G1] = 1;
            // 2

            //--------- o is current node, v is (i2, j2, k2)----------
            //    * * * *
            //   * o v *
            //  * * * *
            //--------------------------------------------------------

            i2 = i + 1;
            j2 = j;
            it = i2;
            jt = j2;
            Boundary(i2, j2);
            if (i2 != it) { BoundaryX[2][G1] = true; }
            else { BoundaryX[2][G1] = false; }
            if (j2 != jt) { BoundaryY[2][G1] = true; }
            else { BoundaryY[2][G1] = false; }
            Gcon[2][G1] = 1 + j2 * Lx + i2;
//           if (j<=Ly/2-5) Kcon[2][G1]=1e-1;
//           else if (j>Ly/2+5) Kcon[2][G1]=1e1;
//           else{
//               Kcon[2][G1]=1e-1*pow(10, (j-Ly/2+5)*0.2);
//            }
            Kcon[2][G1] = 1;
        }
    }
    return;
}


void CrossLinks(double x[], double xm[], /*double xm2[],*/ double theta/*,
                double theta2*/)   // initialize the coordination of the network
{
    int N = Lx * Ly, G1;
    for (int i = 0; i <= (Lx - 1); i++) {
        for (int j = 0; j <= (Ly - 1); j++) {
            G1 = 1 + j * Lx + i;
            x[G1] = i * r[2][0] + (j) * r[1][0];
            x[G1 + N] = i * r[2][1] + j * r[1][1];
        }
    }

    x[2 * N + 1] = xm[0];
    x[2 * N + 2] = xm[1];
    x[2 * N + 3] = theta;

//    x[2 * N + 4] = xm2[0];
//    x[2 * N + 5] = xm2[1];
//    x[2 * N + 6] = theta2;


    return;

}


double MeasureZ(int **Gcon)   //a function measuring network connectivity z//
{
    int i, j, con;

    double z = 0;

    for (i = 0; i <= (Lx - 1); i++) {
        for (j = 0; j <= (Ly - 1); j++) {
            for (con = 0; con <= 2; con++) {
                if (Gcon[con][1 + j * Lx + i] != -1)
                    z = z + 2.0;
            }
        }
    }

    return z / Lx / Ly;
}

void Dilute(int **Gcon, double z) {
    int i, j, con;

    if (MeasureZ(Gcon) < z) cout << "Dilution fails, try a smaller p value!" << endl;

    while (MeasureZ(Gcon) > z) {
        i = rand() % Lx;
        j = rand() % Ly;
        con = rand() % 3;
        Gcon[con][1 + j * Lx + i] = -1;
    }

    return;
}

int ad_direction(double x[], double xm[], double theta) {   //direc tells if any part of cell get out of the net
    int N = Lx * Ly, G1;
    int direc = 0;

    for (int i = 0; i <= (Lx - 1); i++) {
        for (int j = 0; j <= (Ly - 1); j++) {
            double rx, ry;

            G1 = 1 + j * Lx + i;

            rx = x[G1];
            ry = x[G1 + N];

            if (ifcell(rx, ry, xm[0], xm[1], a, b, theta)) {
                if (i == 0 || i == (Lx - 1) || j == 0 || j == (Ly - 1)) {
                    direc = 1;
                }
            }

        }
    }
    return direc;
}

int adhesion_site(double x[], double xm[], double theta) {   //form adhesion site on the lattice site

    int N = Lx * Ly, G1;

    int num = totalnum;   //num=-1 in no adhesion case
    int direc;

    direc = ad_direction(x, xm, theta);  //calculate direc[2]

    cout << "direc= " << direc << endl;

    for (int i = 0; i <= (Lx - 1); i++) {
        for (int j = 0; j <= (Ly - 1); j++) {
            double rx, ry;
            G1 = 1 + j * Lx + i;

            rx = x[G1];
            ry = x[G1 + N];

            int istrue = 0;

            for (int k = 1;
                 k <= totalnum; k++) {   //search if the current node is already connected with an adhesion site
                if (AdhesionIndex[k][0] == G1) {
                    istrue = 1;
                    break;
                }
            }

            double shift[2];
            shift[0] = 0.0;
            shift[1] = 0.0;


            if (ifcell_periodic(rx, ry, xm[0], xm[1], a, b, theta, shift, direc) && (istrue == 0) &&
                (fRand(0, 1) < onrate)) {

                double k_x, p;

                rx = rx - xm[0] + shift[0];
                ry = ry - xm[1] + shift[1];   //this is the coordinate in the ground reference

//                cout<<"rx= "<<rx<<endl;
//                cout<<"ry= "<<ry<<endl;

                double rx2, ry2;        //this is the coordinate in the elliptical reference
                rx2 = rx * cos(theta) + ry * sin(theta);
                ry2 = rx * (-sin(theta)) + ry * cos(theta);


                k_x = k_b - (k_b - k_f) * (rx2 + a) / (2 * a);

                //now in the begining of the contraction cycle, all possible FAs will be formed
                Adhesion[num + 1][0] = rx;
                Adhesion[num + 1][1] = ry;
                Adhesion_shift[num + 1][0] = shift[0];
                Adhesion_shift[num + 1][1] = shift[1];
                AdhesionIndex[num + 1][0] = G1;
                AdhesionIndex[num + 1][1] = G1 + N;
                num = num + 1;

            }
        }
    }

    int totalnum = num;

    for (num = 1; num <= totalnum; num++) {
        Adhesion0[num][0] = Adhesion[num][0] + xm[0];   //Adhesion0 is the original position of adhesion site
        Adhesion0[num][1] = Adhesion[num][1] + xm[1];
        ifremove[num] = 0;                            //all sites remain attached in the beginning
    }

    return totalnum;
}

/*
int adhesion_site2(double x[], double xm[], double theta) {   //form adhesion site on the lattice site

    int N = Lx * Ly, G1;

    int num = totalnum2;   //num=-1 in no adhesion case
    int direc;

    direc = ad_direction(x, xm, theta);  //calculate direc[2]

    cout << "direc= " << direc << endl;

    for (int i = 0; i <= (Lx - 1); i++) {
        for (int j = 0; j <= (Ly - 1); j++) {
            double rx, ry;
            G1 = 1 + j * Lx + i;

            rx = x[G1];
            ry = x[G1 + N];

            int istrue = 0;

            for (int k = 1;
                 k <= totalnum2; k++) {   //search if the current node is already connected with an adhesion site
                if (AdhesionIndex2[k][0] == G1) {
                    istrue = 1;
                    break;
                }
            }

            double shift[2];
            shift[0] = 0.0;
            shift[1] = 0.0;


            if (ifcell_periodic(rx, ry, xm[0], xm[1], a, b, theta, shift, direc) && (istrue == 0) &&
                (fRand(0, 1) < onrate)) {

                double k_x, p;

                rx = rx - xm[0] + shift[0];
                ry = ry - xm[1] + shift[1];   //this is the coordinate in the ground reference

//                cout<<"rx= "<<rx<<endl;
//                cout<<"ry= "<<ry<<endl;

                double rx2, ry2;        //this is the coordinate in the elliptical reference
                rx2 = rx * cos(theta) + ry * sin(theta);
                ry2 = rx * (-sin(theta)) + ry * cos(theta);


                k_x = k_b - (k_b - k_f) * (rx2 + a) / (2 * a);

                //now in the begining of the contraction cycle, all possible FAs will be formed
                Adhesion2[num + 1][0] = rx;
                Adhesion2[num + 1][1] = ry;
                Adhesion_shift2[num + 1][0] = shift[0];
                Adhesion_shift2[num + 1][1] = shift[1];
                AdhesionIndex2[num + 1][0] = G1;
                AdhesionIndex2[num + 1][1] = G1 + N;
                num = num + 1;

            }
        }
    }

    int totalnum2 = num;

    for (num = 1; num <= totalnum2; num++) {
        Adhesion02[num][0] = Adhesion2[num][0] + xm[0];   //Adhesion0 is the original position of adhesion site
        Adhesion02[num][1] = Adhesion2[num][1] + xm[1];
        ifremove2[num] = 0;                            //all sites remain attached in the beginning
    }

    return totalnum2;
}
*/


bool isinNet(double x[], double xm[], double theta) {
    //calculate if the cell has moved out of the network

    int N = Lx * Ly, G1;
    bool flag = false;


    for (int i = 0; i <= (Lx - 1); i++) {
        for (int j = 0; j <= (Ly - 1); j++) {
            double rx, ry;
            G1 = 1 + j * Lx + i;

            rx = x[G1];
            ry = x[G1 + N];

            if (ifcell(rx, ry, xm[0], xm[1], a, b, theta)) {
                flag = true;
                return flag;
            }

        }
    }
    return flag;
}

void shiftcell(double x[], double xm[], double theta) {

    if (isinNet(x, xm, theta)) {
        return;
    } else {

        if (totalnum == 0) {   //no focal adhesion, but the cell is out of the network

            double xc = xm[0];
            double yc = xm[1];

            if (yc > h) {
                yc -= h;
                xc -= Ly / 2;

                double tempX = xc - Ly / 2 * (yc / h);     //if tempX belongs [0, Lx], it is within the network

                if (tempX < 0) {
                    xc += Lx;
                } else if (tempX > Lx) {
                    xc -= Lx;
                }
            } else if (yc < 0) {
                yc += h;
                xc += Ly / 2;

                double tempX = xc - Ly / 2 * (yc / h);

                if (tempX < 0) {
                    xc += Lx;
                } else if (tempX > Lx) {
                    xc -= Lx;
                }
            } else {
                double tempX = xc - Ly / 2 * (yc / h);     //if tempX belongs [0, Lx], it is within the network

                if (tempX < 0) {
                    xc += Lx;
                } else if (tempX > Lx) {
                    xc -= Lx;
                }
            }

            xm[0] = xc;
            xm[1] = yc;

            return;
        }


        for (int i = 1; i <= totalnum; i++) {
            double shiftx = Adhesion_shift[i][0];
            double shifty = Adhesion_shift[i][1];

            double xm0[2];
            xm0[0] = xm[0];
            xm0[1] = xm[1];

            if (fabs(shiftx) > 1e-2 || fabs(shifty) > 1e-2) {
                xm0[0] = xm0[0] - shiftx;
                xm0[1] = xm0[1] - shifty;

                if (isinNet(x, xm0, theta)) {
                    xm[0] = xm0[0];
                    xm[1] = xm0[1];

                    for (int j = 1; j <= totalnum; j++) {
                        Adhesion_shift[j][0] = Adhesion_shift[j][0] - shiftx;
                        Adhesion_shift[j][1] = Adhesion_shift[j][1] - shifty;
                    }

                    return;
                }
            }


        }


    }


}


/*
void shiftcell2(double x[], double xm2[], double theta2) {

    if (isinNet(x, xm2, theta2)) {
        return;
    } else {

        if (totalnum2 == 0) {   //no focal adhesion, but the cell is out of the network

            double xc = xm2[0];
            double yc = xm2[1];

            if (yc > h) {
                yc -= h;
                xc -= Ly / 2;

                double tempX = xc - Ly / 2 * (yc / h);     //if tempX belongs [0, Lx], it is within the network

                if (tempX < 0) {
                    xc += Lx;
                } else if (tempX > Lx) {
                    xc -= Lx;
                }
            } else if (yc < 0) {
                yc += h;
                xc += Ly / 2;

                double tempX = xc - Ly / 2 * (yc / h);

                if (tempX < 0) {
                    xc += Lx;
                } else if (tempX > Lx) {
                    xc -= Lx;
                }
            } else {
                double tempX = xc - Ly / 2 * (yc / h);     //if tempX belongs [0, Lx], it is within the network

                if (tempX < 0) {
                    xc += Lx;
                } else if (tempX > Lx) {
                    xc -= Lx;
                }
            }

            xm2[0] = xc;
            xm2[1] = yc;

            return;
        }


        for (int i = 1; i <= totalnum2; i++) {
            double shiftx = Adhesion_shift2[i][0];
            double shifty = Adhesion_shift2[i][1];

            double xm0[2];
            xm0[0] = xm2[0];
            xm0[1] = xm2[1];

            if (fabs(shiftx) > 1e-2 || fabs(shifty) > 1e-2) {
                xm0[0] = xm0[0] - shiftx;
                xm0[1] = xm0[1] - shifty;

                if (isinNet(x, xm0, theta2)) {
                    xm2[0] = xm0[0];
                    xm2[1] = xm0[1];

                    for (int j = 1; j <= totalnum2; j++) {
                        Adhesion_shift2[j][0] = Adhesion_shift2[j][0] - shiftx;
                        Adhesion_shift2[j][1] = Adhesion_shift2[j][1] - shifty;
                    }

                    return;
                }
            }


        }


    }


}
*/
