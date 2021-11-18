#ifndef ADHESION
#define ADHESION

#include "global_defs.hpp"
#include "params.hpp"
#include "utils.hpp"
#include "geom.hpp"
#include "cell.hpp"
#include "network.hpp"

// Create Random adhesion sites
void Adhesion_sites(Cells &cell_sys, Network &net) {

  // temporary variables
  double x, y;
  int ind = 0;

  // Iterate through all the cells
  for (auto &cell : cell_sys.cells) {

    // Randomly spawn 'Nadh'# of points inside a rectangle
    // of lengths a and b with it's center at the origin
    ind = 0;
    for (int i = 0; i < cell.Nadh; i++) {
      x = randDouble(-1*cell.a/2, cell.a/2);
      y = randDouble(-1*cell.b/2, cell.b/2);

      // Check if the point is inside the ellipse
      if (pow(x/a, 2) + pow(y/b, 2) < 1) {

        // Transform the point so that it is inside the
        // current cell
        cell.Adhesion(ind, 0) = cell.x + x*cos(cell.theta)
                                -y*sin(cell.theta);
        cell.Adhesion(ind, 1) = cell.y + x*sin(cell.theta)
                                +y*cos(cell.theta);
        cell.Adhesion0(ind, 0) = cell.Adhesion(ind, 0);
        cell.Adhesion0(ind, 1) = cell.Adhesion(ind, 1);

        ind++;
      }
    }
  }
}

#endif
