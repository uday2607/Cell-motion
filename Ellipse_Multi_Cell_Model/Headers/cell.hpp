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
    if ((p1[0] < 0) && (p2[0] < 0) &&
        (p3[0] < 0) && (p4[0] < 0)) {

        // push x-coordinate by a
        x = x + a;
    }

    if ((p1[0] > Lx) && (p2[0] > Lx) &&
        (p3[0] > Lx) && (p4[0] > Lx))  {

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
    int Nol, Nadh, index;
    double x, y, theta;
    darray Adhesion, Adhesion0;
    ivect Overlap;
    bool Contract;
    Cell() {};
};

/* Class for Cells in the system */
class Cells {
  public:
    int Ncells;
    int a, b, Nol, Nadh, Lx, Ly;
    std::vector<Cell> cells;
    Cells(int, int, int, double,
          double, int, int);
    bool noCollision(std::vector<Cell>, double, double,
                    double, double, double);
    Contraction(double, double);
};

// Check if two ellipses are intersecting
bool Cells::noCollision(std::vector<Cell> cells, double a, double b,
             double x, double y, double theta) {

    double x1, y1, theta1, a1, b1;

    for (auto cell : cells) {

      // Ellipse 2
      x1, y1 = cell.x, cell.y;
      theta1 = cell.theta;
      a1, b1 = cell.a, cell.b;

      // If they intersect
      if (ellipse_intersects(a, b, x, y, theta,
                             a1, b1, x1, y1, theta1)) {
        return false;
      }
    }
    return true;
}

// Cells class initializer
Cells::Cells(int ncells, int lx, int ly,
             double a, double b, int nol, int N) {

  Ncells = ncells;
  Lx = lx;
  Ly = ly;
  a = a;
  b = b;
  Nol = nol;
  Nadh = N;

  // temp variables
  int i = 0;
  double x, y, theta;

  while (i < Ncells) {

    x = randDouble(0, Lx);
    y = randDouble(0, Ly);
    theta = randDouble(0, 2.*M_PI);

    // Push the cell in to the network
    intoNet(Lx, Ly, a, b, x, y, theta);

    // Check if the Ellipse is clashing with any existing
    // ellipse
    if (noCollision(cells, a, b, x, y, theta)) {
      cells.push_back(Cell());
      cells[i].a, cells[i].b = a, b;
      cells[i].Nol = nol;
      cells[i].Nadh = N;
      cells[i].index = i+1;
      cells[i].x, cells[i].y = x, y;
      cells[i].theta = theta;
      cells[i].Adhesion = -1*darray::Ones(N, 2);
      cells[i].Adhesion0 = -1*darray::Ones(N, 2);
      cells[i].Overlap = -1*ivect::Ones(N);
      cells[i].Contract = true;
      i = i + 1;
    }
  }
}

// Check if a cell overlaps with any other cell
void CheckOverlap(Cell &cell, std::vector<Cell> &cells) {

  // temporary variables
  int NUM = 1000, index, ind;
  double x, y, theta, a, b;
  double x1, y1, theta1, a1, b1;
  double xi, yi;

  dvect phi(NUM);
  darray ell1(NUM, 2);
  darray ell2(NUM, 2);

  // Ellipse 1
  x, y = cell.x, cell.y;
  theta = cell.theta;
  a, b = cell.a, cell.b;
  index = cell.index;

  for (int i = 1; i <= cells.size();  i++) {

    if (i == index+1) {
      continue;
    }

    // Ellipse 2
    x1, y1 = cells[i].x, cells[i].y;
    theta1 = cells[i].theta;
    a1, b1 = cells[i].a, cells[i].b;

    // If they intersect
    if (ellipse_intersects(a, b, x, y, theta,
                           a1, b1, x1, y1, theta1)) {
      ind = Search_ivect_int(cells[i].Overlap, -1);
      cells[i].Overlap[ind] = index;
      ind = Search_ivect_int(cell.Overlap, -1);
      cell.Overlap[ind] = index;
    }
  }
}

// Cell Contraction
void Cells::Contraction(double dt, double step) {

  // temporary variables
  double x, y, c, E;

  // Contract every cell
  for (auto cell : cells) {
    // Check if the cell is overlapping with any
    // other cell
    CheckOverlap(cell, cells);
    if ((cell.Overlap == -1).all()) {
      // Cell doesn't overlap
      // Change semi major axes
      cell.a = cell.a*(1 - lambda*dt/tau);

      // Shift Focal Adhesions normally
      c = lambda*dt/tau;
      for (int i = 0; i < Nadh; i++) {
        // Check if Adhesion exists
        if (cell.Adhesion(i,0) != -1) {
          x = cell.Adhesion(i,0);
          y = cell.Adhesion(i,1);

          cell.Adhesion(i,0) = x - c*pow(cos(cell.theta), 2)*x
                                    - c*cos(cell.theta)
                                       *sin(cell.theta)*y;
          cell.Adhesion(i,1) = y - c*pow(sin(cell.theta), 2)*y
                                    - c*cos(cell.theta)
                                       *sin(cell.theta)*x;
        }
      }

      // Check if the cell semi major axes reached it's min
      // value
      if (cell.a <= lambda*a) {
        cell.Contract = false;
      }
    } else {

      // The cell is overlapping with other cell
      // Calculate the tension created by other cells
      // Along the major axes direction


    }
  }

}

#endif
