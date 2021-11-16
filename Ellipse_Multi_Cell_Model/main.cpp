// Custom Libraries
#include "Headers/global_defs.hpp"
#include "Headers/params.hpp"
#include "Headers/cell.hpp"
#include "Headers/network.hpp"

int main(int argc, char const *argv[]) {

  Lx = 30;
  Ly = 30;
  a = 3;
  b = 2;
  int Ncells = 10;
  int nol = 32;
  int N = 32;

  // Important Variables
  nvariables = 1 + 2*Lx*Ly + 3*Ncells;
  dvect x(nvariables);

  // Initialize the network
  Network net (Lx, Ly);
  net.BuildNetwork();

  // Intialize the Cell System
  Cells cell_sys (Ncells, Lx, Ly, a, b,
                  nol, N);

  for (auto cell : cell_sys.cells) {
    std::cout << cell.index << "\n";
  }

  for (int i = 0; i < Ncells;  i++) {
    delete &cell_sys.cells[i];
  }
  cell_sys.cells.clear();

  return 0;
}
