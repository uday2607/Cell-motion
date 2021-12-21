#ifndef GEOM
#define GEOM

#include "global_defs.hpp"

// Check if a point is in the ellipse
bool inellipse(double a, double b, double x, double y,
               double theta, double xi, double yi) {

    double dis = 0;
    dis = ((pow(((x-xi)*cos(theta)+(y-yi)*sin(theta)),2)/a*a)
          + pow(((x-xi)*sin(theta)-(y-yi)*cos(theta)),2)/b*b);

    return dis <= 1;
}

// Find intersection of two lines
dvect LL_intersection(dvect point1, dvect point2, dvect point3, dvect point4) {

    // return point
    dvect point(2);

    // Line AB represented as a1x + b1y = c1
    double a1 = point2(1) - point1(1);
    double b1 = point1(0) - point2(0);
    double c1 = a1*(point1(0)) + b1*(point1(1));

    // Line CD represented as a2x + b2y = c2
    double a2 = point4(1) - point3(1);
    double b2 = point3(0) - point4(0);
    double c2 = a2*(point3(0))+ b2*(point3(1));

    double determinant = a1*b2 - a2*b1;

    if (determinant == 0)
    {
        // The lines are parallel. This is simplified
        // by returning a pair of LONG_MAX
        point(0) = LONG_MAX;
        point(1) = LONG_MAX;
        return point;
    } else {
        point(0) = (b2*c1 - b1*c2)/determinant;
        point(1) = (a1*c2 - a2*c1)/determinant;
        return point;
    }
}

// Check if two ellipses intersect
bool ellipse_intersects(double a, double b, double x, double y, double theta,
                  double a1, double b1, double x1, double y1, double theta1) {

    // temporary variables
    int NUM = 1000;

    dvect phi(NUM);
    for (int i = 0; i < NUM; i++) {
      phi(i) = i;
    }

    phi = 2.0*M_PI*phi/NUM;
    darray ell1(NUM, 2);
    darray ell2(NUM, 2);

    // Ellipse 1
    ell1(Eigen::seqN(0, NUM), 0) = x + a*cos(theta)*phi.cos()
                                  - b*sin(theta)*phi.sin();
    ell1(Eigen::seqN(0, NUM), 1) = y + a*cos(theta)*phi.sin()
                                  - b*sin(theta)*phi.cos();

    double xi, yi;

    for (int i = 0; i < NUM; i++) {
      xi = ell1(i, 0);
      yi = ell1(i, 1);

      if (inellipse(a1, b1, x1, y1, theta1, xi, yi)) {
        return true;
      }
    }
    return false;
}

#endif
