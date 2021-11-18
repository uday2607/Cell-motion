#ifndef NETWORK_H
#define NETWORK_H

#include "global_defs.hpp"

// Network Class
class Network {
  public:
    int Lx, Ly;
    darray Adh_sites;
    Network (int, int);
};

// Initializer
Network::Network(int lx, int ly) {
  Lx = lx;
  Ly = ly;
}

#endif
