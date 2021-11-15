// Custom Libraries
#include "Headers/global_defs.hpp"
#include "params.hpp"
#include "Headers/cell.hpp"
#include "Headers/network.hpp"

int main(int argc, char const *argv[]) {

  int Ncells = 2;

  // Important Variables
  nvariables = 1 + 2*Lx*Ly + 3*Ncells;
  dvect x(nvariables);

  // Initialize the network
  Network net (Lx, Ly);
  net.BuildNetwork();

  // Intialize the Cell System

  return 0;
}
