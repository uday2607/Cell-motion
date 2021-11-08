// Network functions
#include <cmath>
#include "globals.hpp"
#include "ifcell.hpp"

// Positive modulus
inline int positive_modulo(int i, int n) {
    return (i % n + n) % n;
}

// Impose boundary conditions
void boundary(int &x, int &y) {

  x = positive_modulo(x, Lx);
  y = positive_modulo(y, Ly);

}

// Build network connections
void BuildConnections(iarray &Gcon, darray &Kcon,
                      barray &BoundaryX, barray &BoundaryY) {

  int i2, j2;
  int it, jt;
  int G1;

  // Bulk Connection
  for (int i = 0; i < Lx; i++) {
    for (int j = 0; j < Ly; j++) {
      // In the triangular lattice:
      // the index of each node is 1+Lx*y+x;
      G1 = 1 + j*Lx + i;

      // Store if the neighboring cells
      // are actually in the lattice

      // IN THE XY PLAIN
      //--------- o is current node, v is (i2, j2)----------
      //    v * * *
      //   * o * *
      //  * * * *
      //--------------------------------------------------------
      i2 = i - 1;
      j2 = j + 1;
      it = i2;
      jt = j2;
      boundary(i2, j2);

      //Check if they are the same
      if (i2 != it) {
        BoundaryX(0, G1) = true;
      } else {
        BoundaryX(0, G1) = false;
      }
      if (j2 != jt) {
        BoundaryY(0, G1) = true;
      } else {
        BoundaryY(0, G1) = false;
      }

      //Store the lattice point
      Gcon(0, G1) = 1 + Lx*j2 + i2;
      Kcon(0, G1) = 1;

      //--------- o is current node, v is (i2, j2)----------
      //    * v * *
      //   * o * *
      //  * * * *
      //--------------------------------------------------------
      i2 = i;
      j2 = j + 1;
      it = i2;
      jt = j2;
      boundary(i2, j2);

      //Check if they are the same
      if (i2 != it) {
        BoundaryX(1, G1) = true;
      } else {
        BoundaryX(1, G1) = false;
      }
      if (j2 != jt) {
        BoundaryY(1, G1) = true;
      } else {
        BoundaryY(1, G1) = false;
      }

      //Store the lattice point
      Gcon(1, G1) = 1 + Lx*j2 + i2;
      Kcon(1, G1) = 1;

      //--------- o is current node, v is (i2, j2)----------
      //    * * * *
      //   * o v *
      //  * * * *
      //--------------------------------------------------------
      i2 = i + 1;
      j2 = j;
      it = i2;
      jt = j2;
      boundary(i2, j2);

      //Check if they are the same
      if (i2 != it) {
        BoundaryX(2, G1) = true;
      } else {
        BoundaryX(2, G1) = false;
      }
      if (j2 != jt) {
        BoundaryY(2, G1) = true;
      } else {
        BoundaryY(2, G1) = false;
      }

      //Store the lattice point
      Gcon(2, G1) = 1 + Lx*j2 + i2;
      Kcon(2, G1) = 1;
    }
  }
}

// Function to measure the network connectivity
double computeZ(iarray &Gcon) {

  double z = 0;

  for (int i = 0; i < Lx; i++){
    for (int j = 0; j < Ly; j ++) {
      for (int con = 0; con <= 2; con++) {
        if (Gcon(con, 1+j*Lx+i) != -1) {
          z += 2.0;
        }
      }
    }
  }
  return z/(Lx*Ly);
}

// Function to dilute the network
void NetworkDilute(iarray &Gcon, double z) {

  int i, j, con;

  if (computeZ(Gcon) < z) {
    std::cout << "Dilution failed. Try for a lower value of p" << std::endl;
    exit(2);
  }

  while (computeZ(Gcon) > z) {
    i = randInt(0, Lx);
    j = randInt(0, Ly);
    con = randInt(0, 3);
    // Detach the bond
    Gcon(con, 1+j*Lx+i)  =-1;
  }
}

// Function to initialize the coordination of the network and cell
void CrossLinks(dvect &x, const darray &r, const dvect &xm, const double theta) {

  int G1, N;
  N = Lx*Ly;

  for (int i = 0; i < Lx; i++) {
    for (int j = 0; j < Ly; j++) {
      G1 = 1 + j*Lx + i;
      x(G1) = i*r(2, 0) + j*r(1, 0);
      x(G1 + N) = i*r(2, 1) + j*r(1, 1);
    }
  }

  x(2*N + 1) = xm(0);
  x(2*N + 2) = xm(1);
  x(2*N + 3) = theta;
}

// Perform Affine and Shear Transformation
// on the lattice network
void AffineShear(dvect &x, const double ShearStrain) {

  int G1, N;
  N = Lx*Ly;

  for (int i = 0; i <= (Lx - 1); i++) {
      for (int j = 0; j <= (Ly - 1); j++) {
          G1 = 1 + j * Lx + i;

          // y = y0 * (1+eps)
          x(G1 + N) = (1 + eps)*x(G1 + N);           // y coordinates


          // x = x0 + y*strain
          x(G1) = x(G1) + x(G1 + N)*ShearStrain;    // x coordinates
        }
    }
}

// Function to tell if a cell goes out
// of the network
int adh_direc(const dvect &x, const dvect &xm, const double theta) {

  int N, G1, direc;
  double rx, ry;

  N = Lx*Ly;
  direc = 0;

  for (int i = 0; i < Lx; i++) {
    for (int j = 0; j < Ly; j++) {

      // vortex index
      G1 = 1 + j*Lx + i;

      // vortex coordinates
      rx = x(G1);
      ry = x(G1 + N);

      if (if_node_in_cell(rx, ry, xm(0), xm(1), a, b, theta)) {
        if (i == 0 || i == (Lx - 1) || j == 0 || j == (Ly - 1)) {
            direc = 1;
        }
      }
    }
  }
  return direc;
}

// Form a FA on the Cell from lattice site
int adhesion_site(const dvect &x, const dvect &xm, const double theta) {

  int N, G1, num, direc;
  N = Lx*Ly;
  num = totalnum; // num = -1 in no adhesion case

  // Calculate the direction
  direc = adh_direc(x, xm, theta);

  std::cout << " Direc = " << direc << std::endl;

  for (int i = 0; i < Lx; i++) {
      for (int j = 0; j < Ly; j++) {

          double rx, ry;
          // vortex index
          G1 = 1 + j*Lx + i;

          // vortex coordinates
          rx = x(G1);
          ry = x(G1 + N);

          int isconnected = 0;

          // search if the current node is already
          // connected with an adhesion site
          for (int k = 1; k <= totalnum; k++) {
              if (AdhesionIndex(k, 0) == G1) {
                  isconnected = 1;
                  break;
              }
          }

          // temporary variable
          double shift[2];
          shift[0] = 0.0;
          shift[1] = 0.0;

          if (if_cell_is_periodic(rx, ry, xm(0), xm(1), a, b,
                                  theta, shift, direc)
              && (isconnected == 0) && (randDouble(0, 1) < onrate)) {

              double k_x;

              // this is the coordinate in the lab frame
              rx = rx - xm(0) + shift[0];
              ry = ry - xm(1) + shift[1];

              // this is the coordinate in the elliptical reference
              double rx2, ry2;
              rx2 = rx*cos(theta) + ry*sin(theta);
              ry2 = rx*(-sin(theta)) + ry*cos(theta);

              // Front attachments > back
              k_x = k_b - (k_b - k_f)*(rx2 + a)/(2*a);

              // Now in the begining of the contraction cycle,
              // all possible FAs will be formed
              Adhesion(num + 1, 0) = rx;
              Adhesion(num + 1, 1) = ry;
              Adhesion_shift(num + 1, 0) = shift[0];
              Adhesion_shift(num + 1, 1) = shift[1];
              AdhesionIndex(num + 1, 0) = G1;
              AdhesionIndex(num + 1, 1) = G1 + N;
              num = num + 1;
            }
          }
        }

  // Total number of adhesion sites
  int totalnum = num;

  for (num = 1; num <= totalnum; num++) {
    // Adhesion0 is the original position of adhesion site
    Adhesion0(num ,0) = Adhesion(num ,0) + xm(0);
    Adhesion0(num ,1) = Adhesion(num ,1) + xm(1);
    // all sites remain attached in the beginning
    ifremove(num) = 0;
  }

  return totalnum;
}
