#ifndef CELL
#define CELL

#include "global_defs.hpp"
#include "params.hpp"
#include "utils.hpp"
#include "geom.hpp"

class Cell {
    int a, b, nol, N, index;
    double x, y, theta;
    darray Adhesion, Adhesion0;

  public:
    Cell(int, int, int, int, int, double, double, double);
};

Cell::Cell(int a, int b, int nol, int N, int index,
           double x, double y, double theta) {

  a = a;
  b = b;
  nol = nol;
  N = N;
  index = index;
  x = x;
  y = y;
  theta = theta;
}

class Cells {
    int Ncells;
    int a, b, nol, N, Lx, Ly;
    std::vector<Cell> cells;
  public:
    Cells(int);
};

Cells::Cells(int ncells, int lx, int ly,
             int a, int b, int nol, int N) {

  Ncells = ncells;
  Lx = lx;
  Ly = ly;
  a = a;
  b = b;
  nol = nol;
  N = N;

  // temp variables
  double x, y, theta;

  for (int i = 0; i < Ncells; i++) {

    // Initial Cell
    if (i == 0) {
      x = randDouble(0, Lx);
      y = randDouble(0, Ly);

      // Shear Transform square to rhombus
      // **** ->    ****
      // **** ->   ****
      // **** ->  ****
      // **** -> ****

      theta = randDouble(0, 2.*M_PI);

      // Check if the cell is in the network

      cells.push_back(Cell(a, b, nol, N, i+1,
                          x, y, theta));

    }
  }

}

#endif
