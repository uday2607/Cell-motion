#include <cmath>
#include "globals.hpp"

// Contraction in Cell
// Breakage of FAs

void contraction(darray &Adhesion, double theta, int fold) {

  for (int i = 0; i <= N; i++) {
    double x, y;
    double c0 = lambda*dt*((double) fold)/tau;

    // Previous FA
    x = Adhesion(i, 0);
    y = Adhesion(i, 1);

    // Updated FA
    Adhesion(i, 0) = x - c0*cos(theta)*cos(theta)*x - c0*sin(theta)*cos(theta)*y;
    Adhesion(i, 1) = y - c0*sin(theta)*cos(theta)*x - c0*sin(theta)*sin(theta)*y;
  }
}
