#include "utils.hpp"
#include "mnbrak.hpp"
//Tolerance passed to dbrent.
#define TOL 2.0e-8

void dlinmin(double p[], double xi[], int n, double *fret,
             double (*func)(double []), void (*dfunc)(double [], double [])) {

//Given an n-dimensional point p[1..n] and an n-dimensional direction xi[1..n], moves and
//resets p to where the function func(p) takes on a minimum along the direction xi from p,
//and replaces xi by the actual vector displacement that p was moved. Also returns as fret
//the value of func at the returned location p. This is actually all accomplished by calling the
//routines mnbrak and dbrent.

  double xx, xmin, fx, fb, fa, bx, ax;
  int ncom = n;

  dvect pcom(n+1);
  dvect xicom(n+1);

  for (int j = 1; j <= n; j++){
    pcom(j) = p(j);
    xicom(j) = xi(j);
  }

  // Initial guess for brackets
  ax = 0.0;
  xx = 1.0;

  mnbrak(&ax, &xx, &bx, &fa, &fx, &fb, f1dim);
  *fret = dbrent(ax, xx, bx, f1dim, dfidim)

}
