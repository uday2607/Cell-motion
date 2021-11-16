#ifndef NETWORK_H
#define NETWORK_H

#include "global_defs.hpp"

// Network Class
class Network {
  public:
    int Lx, Ly;
    darray Coords;
    iarray Connections;
    Network (int, int);
    void BuildNetwork();
};

Network::Network(int lx, int ly) {
  Lx = lx;
  Ly = ly;

  // Initialize the arrays
  Coords = darray(Lx*Ly, 2);
  Connections = iarray::Zero(Lx*Ly, 2);
}

void Network::BuildNetwork() {

  int ind;
  double sqrt_inv_three = sqrt(1/3);

  // Bulk Connection
  for (int i = 0; i < Lx; i++) {
    for (int j = 0; j < Ly; j++) {
      ind = i + Lx*j;

      // Triangular lattice
      Coords(ind, 0) = (2*i+j)*sqrt_inv_three;
      Coords(ind, 1) = j;
    }
  }
}

#endif
