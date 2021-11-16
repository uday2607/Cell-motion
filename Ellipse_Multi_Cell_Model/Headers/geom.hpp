#ifndef GEOM
#define GEOM

#include "global_defs.hpp"

// Check if a point is in the ellipse
bool inellipse(double x, double y, double theta, double a,
               double b, double xi, double yi) {

    double dis = 0;
    dis = ((pow(((x-xi)*cos(theta)+(y-yi)*sin(theta)),2)/a*a)
          + pow(((x-xi)*sin(theta)-(y-yi)*cos(theta)),2)/b*b);

    return dis <= 1;
}

#endif
