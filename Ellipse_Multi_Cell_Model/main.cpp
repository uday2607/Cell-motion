// Custom Libraries
#include "Headers/global_defs.hpp"
#include "Headers/params.hpp"
#include "Headers/cell.hpp"
#include "Headers/network.hpp"
#include "Headers/adhesion.hpp"

int main(int argc, char const *argv[]) {

  Lx = 30;
  Ly = 30;
  a = 3;
  b = 2;
  int Ncells = 10;
  int nol = 32;
  int N = 4*(floor(a)+1)*(floor(b)+1);

  // Initialize the network
  Network net (Lx, Ly);

  // Intialize the Cell System
  Cells cell_sys (Ncells, Lx, Ly, a, b,
                  nol, N);

  // Create adhesion sites
  Adhesion_sites(cell_sys, net);

  for (auto cell : cell_sys.cells) {
    std::cout << cell.Adhesion << "\n";
  }

  return 0;
}
