#include <cmath>
#include <iostream>
#include <utils.h>
#include "linmin.h"

// Here ITMAX is the maximum allowed number of iterations,
// while EPS is a small number to rectify the special case
// of converging to exactly zero function value.
#define EPS 1.0e-14
#define ITMAX 100000

// Given a starting point p[1..n], Fletcher-Reeves-Polak-Ribiere
// minimization is performed on a function func, using its
// gradient as calculated by a routine dfunc. The convergence tolerance
// on the function value is input as ftol. Returned quantities are p
// (the location of the minimum), iter (the number of iterations that
// were performed), and fret (the minimum value of the function).
// The routine linmin is called to perform line minimizations.
void frprmn(dvect &p, int n, double ftol, int *iter, double *fret,
		double (*func)(dvect), void (*dfunc)(dvect, dvect)) {

    int gg, gam, fp, dgg;

    dvect g(n+1);
    dvect h(n+1);
    dvect xi(n+1);

    // Initialization
    fp = (*func)(p);
    (*dfunc)(p, xi);

    for (int j = 1; j <= n; j++) {
      g(j) = -xi(j);
      xi(h) = h(j) = g(j);
    }

    for (int its = 1; its <= ITMAX; its++) {
      // Loop over iterations
      *iter = its;

      // Perform Line Minimization
      linmin(p, xi, n, fret, func, dfunc);

      fp = *fret;
      (*dfunc)(p, xi);

      double res = xi(Eigen::seq(1, n)).max()
      dgg = gg = 0.0;

      if (res <= 1e-5) {
        std::cout << "Good enough: " << abs(*fret-fp) << endl;
        return;
      }

      for (int j = 1; j <= n; j++) {
        gg += g(j)*g(j);
        // Polak-Ribiere statement
        dgg += (xi(j) + g(j))*xi(j);
      }

      if (*iter%20 == 0) {
        std::cout << "Number of iterations: " << *iter << endl;
        std::cout << "Gradient magnitude: " << res << "||||" << "current energy " << *fret << std::endl;
      }

      if (dgg == 0.0) {
        // Unlikely. If gradient is exacyly zero
        // we are already done.
        std::cout << "Zero Gradient -> No more worrying" << std::endl;
        return;
      }
      gam = dgg/gg;
      for (int j = 1; j <= n; j++) {
        g(j) = -xi(j);
        xi(j) = h(j) = g(j)+gam*h(j);
      }
    }
    std::cout << "Nothing happening..." << std::endl;
    return;
}
