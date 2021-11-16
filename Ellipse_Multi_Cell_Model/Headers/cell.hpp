#ifndef CELL
#define CELL

#include "global_defs.hpp"
#include "params.hpp"
#include "utils.hpp"
#include "geom.hpp"

// Helper functions
void intoNet(int Lx, int Ly, double a, double b,
             double &x, double &y, double theta) {

    // 2D Coordinates
    double p1[2], p2[2], p3[2], p4[2];

    // Compute extremes of the ellipse
    p1[0] = x - a*cos(theta);
    p1[1] = y - a*sin(theta);
    p2[0] = x + a*cos(theta);
    p2[1] = y + a*sin(theta);
    p3[0] = x - b*sin(theta);
    p3[1] = y + b*cos(theta);
    p4[0] = x + b*sin(theta);
    p4[1] = y - b*cos(theta);

    // Check if the cell is in the lattice
    // by checking if the cells are on the same side
    // of the lattice line segments
    if (((p1[1] > sqrt(3)*p1[0]) &&
        (p2[1] > sqrt(3)*p2[0]) &&
        (p3[1] > sqrt(3)*p3[0]) &&
        (p4[1] > sqrt(3)*p4[0]))) {

        // push x-coordinate by a
        x = x + a;
    }

    if (((p1[1] < sqrt(3)*p1[0] - 2*Lx) &&
        (p2[1] < sqrt(3)*p2[0] - 2*Lx) &&
        (p3[1] < sqrt(3)*p3[0] - 2*Lx) &&
        (p4[1] < sqrt(3)*p4[0] - 2*Lx))) {

        // push x-coordinate by a
        x = x - a;
    }

    if ((p1[1] < 0) && (p2[1] < 0) &&
        (p3[1] < 0) && (p4[1] < 0)) {

        // push y-coordinate by b
        y = y + b;
    }

    if ((p1[1] > Ly) && (p2[1] > Ly) &&
        (p3[1] > Ly) && (p4[1] > Ly)) {

        // push y-coordinate by b
        y = y - b;
    }
}
/* Class for Cell in the System */
class Cell {
  public:
    double a, b;
    int nol, N, index;
    double x, y, theta;
    darray Adhesion, Adhesion0;
    Cell(double a, double b, int nol, int N, int index,
               double x, double y, double theta) {

      a = a;
      b = b;
      nol = nol;
      N = N;
      index = index;
      x = x;
      std::cout << x << std::endl;
      y = y;
      theta = theta;
    }
};

/* Class for Cells in the system */
class Cells {
  public:
    int Ncells;
    int a, b, nol, N, Lx, Ly;
    std::vector<Cell> cells;
    Cells(int, int, int, double,
          double, int, int);
};

// Check if two ellipses are intersecting
bool noCollision(std::vector<Cell> cells, double a, double b,
             double x, double y, double theta) {


    // temporary variables
    int NUM = 10000;
    double x1, y1, theta1, a1, b1;
    double xi, yi;

    dvect phi(NUM);
    darray ell1(NUM, 2);
    darray ell2(NUM, 2);

    // Ellipse 1
    ell1(Eigen::seqN(0, NUM), 0) = x + a*cos(theta)*phi.cos()
                                     - b*sin(theta)*phi.sin();
    ell1(Eigen::seqN(0, NUM), 1) = y + a*cos(theta)*phi.sin()
                                     - b*sin(theta)*phi.cos();

    for (auto cell : cells) {
      x1 = cell.x;
      y1 = cell.y;
      theta1 = cell.theta;
      a1 = cell.a;
      b1 = cell.b;

      ell2(Eigen::seqN(0,NUM),0) = x1 + a1*cos(theta1)*phi.cos()
                                      - b1*sin(theta1)*phi.sin();
      ell2(Eigen::seqN(0,NUM),1) = y1 + a1*cos(theta1)*phi.sin()
                                      - b1*sin(theta1)*phi.cos();

      for (int i = 0; i < NUM; i++){
        xi = ell1(i, 0);
        yi = ell1(i, 1);
        if (inellipse(x1, y1, theta1, a1, b1, xi, yi)) {
          return false;
        }
      }

      for (int i = 0; i < NUM; i++){
        xi = ell2(i, 0);
        yi = ell2(i, 1);
        if (inellipse(x, y, theta, a, b, xi, yi)) {
          return false;
        }
      }
    }

    return true;
}

Cells::Cells(int ncells, int lx, int ly,
             double a, double b, int nol, int N) {

  Ncells = ncells;
  Lx = lx;
  Ly = ly;
  a = a;
  b = b;
  nol = nol;
  N = N;

  // temp variables
  int i = 0;
  double x, y, theta;
  double sqrt_inv_t = 1/sqrt(3);
  std::vector<Cell> temp_cells;

  while (i < Ncells) {

    std::cout << i << "\n";

    x = randDouble(0, Lx);
    y = randDouble(0, Ly);

    // Shear Transform square to rhombus
    // **** ->    ****
    // **** ->   ****
    // **** ->  ****
    // **** -> ****
    x = (2*x+y)*sqrt_inv_t;

    theta = randDouble(0, 2.*M_PI);

    // Push the cell in to the network
    intoNet(Lx, Ly, a, b, x, y, theta);

    // Check if the Ellipse is clashing with any existing
    // ellipse
    if (noCollision(cells, a, b, x, y, theta)) {
      Cell temp_cell = Cell(a, b, nol, N, i+1,
                          x, y, theta);
      temp_cells.push_back(temp_cell);
      cells[i] = temp_cell;
      i = i + 1;
    }
  }
}

#endif
