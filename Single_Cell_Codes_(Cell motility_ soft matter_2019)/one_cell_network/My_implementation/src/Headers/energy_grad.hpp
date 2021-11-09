#include "Headers/globals.hpp"

#define tiny 1e-16

// Function to calculate the adhesion energy
double AdEnergy(dvect &x) {

  int N = Lx*Ly;
  ivect IND = Eigen::seq(1, total)

  // cell coordinates
  double xc_new = x(2*N+1);
  double yc_new = x(2*N+2);
  double theta = x(2*N+3);

  // fix one corner of the network to prevent possible shift
  x(1) = 0;
  x(1+N) = 0;

  // Energy to be calculated
  double E = 0;

  // Get the indices of the adhesion sites
  std::vector<int> ind;
  // ifremove[i] == 1 -> FA is removed
  std::for_each(&ifremove(1), &ifremove(total), isEqualInt(&ind, 1));

  // Get the coordinates
  dvect xx = Adhesion(IND, 0);
  dvect yy = Adhesion(IND, 1);

  ivect ix = AdhesionIndex(IND, 0);
  ivect iy = AdhesionIndex(IND, 1);

  // now x0, y0 is no longer fixed
  // shift x0, y0 due to the periodic conditions
  dvect x0 = x(ix) + Adhesion_shift(IND, 0);
  dvect y0 = x(iy) + Adhesion_shift(IND, 1);

  dvect x_new = cos(theta)*xx - sin(theta)*yy;
  dvect y_new = sin(theta)*xx + cos(theta)*yy;

  dvect fx = -(x_new + xc_new - x0);
  dvect fy = -(y_new + yc_new - y0);

  // Compute the energy
  E = 0.5*ks*((fx.square() + fy.square()).sum())
  return E;
}

void Adgrad(dvect &x, dvect &g) {

  int N = Lx*Ly;
  ivect IND = Eigen::seq(1, total)

  // cell coordinates
  double xc_new = x(2*N+1);
  double yc_new = x(2*N+2);
  double theta = x(2*N+3);

  // fix one corner of the network to prevent possible shift
  x(1) = 0;
  x(1+N) = 0;

  // Get the coordinates
  dvect xx = Adhesion(IND, 0);
  dvect yy = Adhesion(IND, 1);

  ivect ix = AdhesionIndex(IND, 0);
  ivect iy = AdhesionIndex(IND, 1);

  // now x0, y0 is no longer fixed
  // shift x0, y0 due to the periodic conditions
  dvect x0 = x(ix) + Adhesion_shift(IND, 0);
  dvect y0 = x(iy) + Adhesion_shift(IND, 1);

  dvect x_new = cos(theta)*xx - sin(theta)*yy;
  dvect y_new = sin(theta)*xx + cos(theta)*yy;

  dvect fx = -(x_new + xc_new - x0);
  dvect fy = -(y_new + yc_new - y0);

  // Compute torque
  dvect torque = -(xc_new - x0)*y_new + (yc_new - y0)*x_new;

  // Store to g
  g(Adhesion_shift(IND, 0)) = g(Adhesion_shift(IND, 0)) + ks*fx;
  g(Adhesion_shift(IND, 1)) = g(Adhesion_shift(IND, 1)) + ks*fy;

  g(2*N+1) = g(2*N+1) - ks*(fx.sum());
  g(2*N+2) = g(2*N+2) - ks*(fy.sum());
  g(2*N+3) = g(2*N+3) + ks*(torque.sum());

  // fix one corner
  g(1) = 0.0;
  g(1+N) = 0.0;
}

// Compute energy change due to the network stretching
double EnergyStretch(dvect &x) {

  int G1, G2, N = Lx*Ly;
  double Energy = 0, dx, dy, dz, dis;

  //fix one corner of the network to prevent possible shift
  x(1) = 0;
  x(1 + N) = 0;

  for (int i = 0; i < Lx; i++) {
      for (int j = 0; j < Ly; j++) {
          G1 = 1 + j*Lx + i;
          for (int con = 0; con < 3; con++) {
              G2 = Gcon(con, G1);
              if (G2 != -1) {
                  //-----beautiful code for phantomalization---------------

                  dx = x(G2) - x(G1);
                  dy = x(G2 + N) - x(G1 + N);

                  //for con==1,2,3 all the following codes are valid,
                  //since they only change one index at a time

                  if (BoundaryX(con, G1)) {
                      dx = dx + Lx * (con == 0 ? -1.0 : 1.0);
                  }
                  if (BoundaryY(con, G1)) {
                      dy = dy + h;
                      dx = dx + 0.5*Ly;
                      //Lees-edwards boundary condition
                  }
                  dis = sqrt(dx*dx + dy*dy);
                  Energy += Kcon(con, G1) * powMy(dis - 1.0, 2)/2.0;
              }
          }
      }
  }
  return Energy;
}

// Compute the energy gradient due to network stretching
void StretchEnergyGrad(dvect &x, dvect &Grad) {

    double dx, dy, temp;
    int G1, G2, G3, N = Lx*Ly;
    double x1, x2, x3, y1, y2, y3, dis;

    //fix one corner of the network to prevent possible shift
    x(1) = 0;
    x(1 + N) = 0;

    //the last three elements are xc, yc and theta
    Grad = dvect::Zero(nvariables);

    for (int i = 0; i < Lx; i++) {
        for (int j = 0; j < Ly; j++) {
            G1 = 1 + j*Lx + i;
            for (int con = 0; con < 3; con++) {
                G2 = Gcon(con, G1);
                if (G2 != -1) {
                    //-----beautiful code for phantomalization---------------

                    dx = x(G2) - x(G1);
                    dy = x(G2 + N) - x(G1 + N);

                    //for con==1,2,3 all the following codes are valid,
                    //since they only change one index at a time

                    if (BoundaryX(con, G1)) {
                        dx = dx + Lx * (con == 0 ? -1.0 : 1.0);
                    }
                    if (BoundaryY(con, G1)) {
                        dy = dy + h;
                        dx = dx + 0.5*Ly;
                        //Lees-edwards boundary condition
                    }
                    dis = sqrt(dx*dx + dy*dy);

                    temp = Kcon(con, G1)*(dis - 1.0)/(dis + tiny);

                    Grad(G1) -= dx*temp;
                    Grad(G1 + 1*N) -= dy*temp;
                    Grad(G2) += dx*temp;
                    Grad(G2 + 1*N) += dy*temp;

                    G3 = Gcon(con, G2);
                    if (G3 != -1) {
                        x1 = x(G1);
                        y1 = x(G1 + N);
                        x2 = x(G2);
                        y2 = x(G2 + N);
                        x3 = x(G3);
                        y3 = x(G3 + N);

                        if (BoundaryX(con, G1)) {
                            x2 += Lx*(con == 0 ? -1 : 1);
                            x3 += Lx*(con == 0 ? -1 : 1);
                        } else if (BoundaryX(con, G2)) {
                            x3 += Lx*(con == 0 ? -1 : 1);
                        }

                        if (BoundaryY(con, G1)) {
                            y2 += h;
                            y3 += h;
                            x2 += 0.5*Ly;
                            x3 += 0.5*Ly;
                        } else if (BoundaryY(con, G2)) {
                            y3 += h;
                            x3 += 0.5*Ly;
                        }

                        double xji, yji, xjk, yjk;

                        xji = x2 - x1;
                        yji = y2 - y1;
                        xjk = x2 - x3;
                        yjk = y2 - y3;

                        double lji, ljk, dot;

                        dot = xji * xjk + yji * yjk;

                        lji = sqrt(xji*xji + yji*yji) + tiny;
                        ljk = sqrt(xjk*xjk + yjk*yjk) + tiny;

                        double temp = dot/(lji*ljk);

                        if (temp < -1.0) {
                            temp = -1.0;
                        } else if (temp > 1.0) {
                            temp = 1.0;
                        }

                        double angle = acos(temp) - M_PI;
                        double common_factor = kappa * angle * (-1 / (sqrt(1 - temp * temp) + tiny));
                        double c1, c2, c3, c4;

                        c1 = xjk/(lji*ljk) - (xji*dot)/(ljk*lji*lji*lji);
                        c2 = xji/(lji*ljk) - (xjk*dot)/(lji*ljk*ljk*ljk);
                        //change x->y in c1;
                        c3 = yjk/(lji*ljk) - (yji*dot)/(ljk*lji*lji*lji);
                        c4 = yji/(lji*ljk) - (yjk*dot)/(lji*ljk*ljk*ljk);

                        Grad(G1) -= common_factor*c1;         //g[x1]
                        Grad(G2) += common_factor*(c1 + c2);  //g[x2]
                        Grad(G3) -= common_factor*c2;         //g[x3]


                        Grad(G1 + 1*N) -= common_factor*c3;         //g[y1]
                        Grad(G2 + 1*N) += common_factor*(c3 + c4);  //g[y2]
                        Grad(G3 + 1*N) -= common_factor*c4;         //g[y3]

                    }
                }
            }
        }
    }

    // fix one corner
    Grad[1] = 0.0;
    Grad[1 + N] = 0.0;
}

// Compute the angle phi
double cal_phi(double x_i, double y_i, double x_j,
               double y_j, double x_k, double y_k) {

    double temp, y;

    temp = ((x_j - x_i)*(x_j - x_k) + (y_j - y_i)*(y_j - y_k))/
           ((sqrt((x_j - x_i)*(x_j - x_i) + (y_j - y_i)*(y_j - y_i)) + tiny)*
            (sqrt((x_j - x_k)*(x_j - x_k) + (y_j - y_k)*(y_j - y_k)) + tiny));

    if (temp < -1.0) {
        temp = -1.0;
    } else if (temp > 1.0) {
        temp = 1.0;
    }

    // Phi
    y = acos(temp) - M_PI;

    return y;
}

// Calculate the Energy due to the bending
double EnergyBend(dvect &x) {

  double Energy = 0, x1, x2, x3, y1, y2, y3;
  int G1, G2, G3, N = Lx * Ly;

  for (int i = 0; i < Lx; i++) {
      for (int j = 0; j < Ly; j++) {
          G1 = 1 + j*Lx + i;
          for (int con = 0; con < 3; con++) {
              G2 = Gcon(con, G1);
              if (G2 != -1) {
                G3 = Gcon(con, G2);
                if (G3 != -1) {
                    x1 = x(G1);
                    y1 = x(G1 + N);
                    x2 = x(G2);
                    y2 = x(G2 + N);
                    x3 = x(G3);
                    y3 = x(G3 + N);

                    if (BoundaryX(con, G1)) {
                        x2 += Lx*(con == 0 ? -1 : 1);
                        x3 += Lx*(con == 0 ? -1 : 1);
                    } else if (BoundaryX(con, G2)) {
                        x3 += Lx*(con == 0 ? -1 : 1);
                    }

                    if (BoundaryY(con, G1)) {
                        y2 += h;
                        y3 += h;
                        x2 += 0.5 * Ly;
                        x3 += 0.5 * Ly;
                    } else if (BoundaryY(con, G2)) {
                        y3 += h;
                        x3 += 0.5 * Ly;
                    }

                    double angle = cal_phi(x1, y1, x2, y2, x3, y3);
                    Energy = Energy + 0.5*angle*angle;
              }
          }
      }
    }
  }
  return Energy;
}

// Compute total energy
double energy(dvect &x) {

  double Energy;
  Energy = AdEnergy(x) + EnergyStretch(x) + kappa*EnergyBend(x);
  return E;
}

// Compute the total gradient
void grad(dvect &x, dvect &g) {

  StretchEnergyGrad(x, g);
  Adgrad(x, g);
}

// Compute Energy gradient (numerical)
void Numerical_EnergyGrad(dvect &x, dvect &Grad) {

    dvect x0(4) = x;
    double energy1, energy2;

    for (int j = 0; j < 4; j++) {

        x0(j) = x(j) + 1e-4;
        energy1 = energy(x0);

        x0(j) = x(j) - 1e-4;
        energy2 = energy(x0);

        // Gradient
        Grad(j) = (energy1 - energy2)/2e-4;

        x0(j) = x(j);
    }
}
