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
    Cell() {};
};

/* Class for Cells in the system */
class Cells {
  public:
    int Ncells;
    int cur_index = -1;
    int cur_phase = 0;
    int a, b, Nol, Nadh, Lx, Ly;
    std::vector<Cell> cells;
    Cells(int, int, int, double,
          double, int, int);
    bool noCollision(std::vector<Cell>, double, double,
                    double, double, double);
    Contraction(double, double, double, double);
};

// Check if two ellipses are intersecting
bool Cells::noCollision(std::vector<Cell> cells, double a, double b,
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
  double sqrt_inv_t = 1/sqrt(3);

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
      cells[i].Adhesion = -darray::Ones(N, 2);
      cells[i].Adhesion0 = -darray::Ones(N, 2);
      cells[i].Overlap = -ivect::Ones(N);
      i = i + 1;
    }
  }
}

// Check if a cell overlaps with any other cell
void CheckOverlap(Cell &cell, std::vector<cell> &cells) {

  // temporary variables
  int NUM = 10000, index, ind;
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
  ell1(Eigen::seqN(0, NUM), 0) = x + a*cos(theta)*phi.cos()
                                   - b*sin(theta)*phi.sin();
  ell1(Eigen::seqN(0, NUM), 1) = y + a*cos(theta)*phi.sin()
                                   - b*sin(theta)*phi.cos();

  for (int i = 1; i <= cells.size();  i++) {

    if (i == index+1) {
      continue;
    }

    x1 = cells[i].x;
    y1 = cells[i].y;
    theta1 = cells[i].theta;
    a1 = cells[i].a;
    b1 = cells[i].b;

    ell2(Eigen::seqN(0,NUM),0) = x1 + a1*cos(theta1)*phi.cos()
                                    - b1*sin(theta1)*phi.sin();
    ell2(Eigen::seqN(0,NUM),1) = y1 + a1*cos(theta1)*phi.sin()
                                    - b1*sin(theta1)*phi.cos();

    for (int i = 0; i < NUM; i++){
      xi = ell1(i, 0);
      yi = ell1(i, 1);
      if (inellipse(x1, y1, theta1, a1, b1, xi, yi)) {
        ind = Search_ivect_int(cell.Overlap, -1);
        cell.Overlap[ind] = cell[i].index;
        ind = Search_ivect_int(cell[i].Overlap, -1);
        cell[i].Overlap[ind] = index;
      }
    }

    for (int i = 0; i < NUM; i++){
      xi = ell2(i, 0);
      yi = ell2(i, 1);
      if (inellipse(x, y, theta, a, b, xi, yi)) {
        ind = Search_ivect_int(cell.Overlap, -1);
        cell.Overlap[ind] = cell[i].index;
        ind = Search_ivect_int(cell[i].Overlap, -1);
        cell[i].Overlap[ind] = index;
      }
    }
  }
}

// Cell Contraction
void Cells::Contraction(double lambda, double dt, double tau,
                        double step) {

  // temporary variables
  double x, y, c;
  Cell temp_cell;

  // Check if any cell is in contraction phase
  if (cur_index == -1) {
    // Start a new contraction
    cur_phase = 1;
    cur_index = randInt(0, Ncells);
    temp_cell = cells[cur_index];

    c = lambda*dt/tau*cur_phase;

    // Check if the cell is overlapping with any
    // other cell
    CheckOverlap(temp_cell, cells);
    if ((temp_cell.Overlap == -1).all()) {
      // Shift Focal Adhesions
      for (int i = 0; i < Nadh; i++) {
        // Check if Adhesion exists
        if (temp_cell.Adhesion(i,0) != -1) {
          x = temp_cell.Adhesion(i,0);
          y = temp_cell.Adhesion(i,1);

          temp_cell.Adhesion(i,0) = x - c*pow(cos(temp_cell.theta), 2)*x
                                    - c*cos(temp_cell.theta)
                                       *sin(temp_cell.theta)*y;
          temp_cell.Adhesion(i,1) = y - c*pow(sin(temp_cell.theta), 2)*y
                                    - c*cos(temp_cell.theta)
                                       *sin(temp_cell.theta)*x;
        }
      }
    } else {
      
    }

    cells[cur_index] = temp_cell;

  } else {

  }

}

#endif
