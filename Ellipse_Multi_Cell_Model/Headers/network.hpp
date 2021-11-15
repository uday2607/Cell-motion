#ifndef NETWORK_H
#define NETWORK_H

#include "global_defs.hpp"

// gives the corresponding points according to
// periodic boundary condition
void Boundary(int &i2, int &j2) {
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

// Network Class
class Network {
    int Lx, Ly;
    iarray Gcon;
    darray Kcon;
    barray BoundaryX;
    barray BoundaryY;
  public:
    Network (int, int);
    void BuildNetwork();
};

Network::Network(int lx, int ly) {
  Lx = lx;
  Ly = ly;
}

void Network::BuildNetwork() {


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
