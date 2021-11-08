// Standard Libraries
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <vector>
#include <random>
#include <eigen3/Eigen/Core>
#include "random/pcg_random.hpp"

// Custom Libraries
#include "Headers/globals.hpp"
#include "Headers/fileio.hpp"
#include "Headers/network.hpp"
#include "Headers/cell.hpp"

int main(int argc, char const *argv[]) {

  //Input
  std::string out, ext, seed;

  /*std::cout << "Output filename: ";
  std::cin >> out;
  std::cout << "Extension: ";
  std::cin >> ext;
  std::cout << "Seed: ";
  std::cin >> seed;*/
  out = "t";
  ext = "10";
  seed = "0";

  //Output
  out_dir = "Output/" + std::string(out);

  // Store Length of the lattice used
  std::cout << Ly << std::endl;
  OutputLx(ext, seed);

  // Normal extension
  eps = std::stof(ext) * 0.20;
  h = Ly*0.86602540378443864676*(1+eps);

  // Local Variables
  iarray Gcon(3, 1+Lx*Ly);
  darray Kcon(3, 1+Lx*Ly);
  barray BoundaryX(3, 1+Lx*Ly);
  barray BoundaryY(3, 1+Lx*Ly);
  darray r(3, 2);

  // Displacement vector
  // xnew=xold+r[con][0], ynew=yold+r[con][1]
  r(0,0) = -0.5;
  r(0,1) = 0.86602540378443864676;

  r(1,0) = 0.5;
  r(1,1) = 0.86602540378443864676;

  r(2,0) = 1;
  r(2,1) = 0;

  dvect xm(2);               // Center of the cell
  xm[0] = double(Lx)*(3.0/4.0) - 10;
  xm[1] = double(Ly)/2 * sqrt(3)/2;

  // Angle of the cell
  double theta = randDouble(0.0, 2.*M_PI);
  double theta0 = theta;

  // Print the initial position and orientation
  std::cout << "The initial position and orientation";
  std::cout << " of the cell is: " << std::endl;
  std::cout << " (x0, y0) = (" << xm[0] << ", " << xm[1] << ")\n";
  std::cout << "theta0 = " << theta << std::endl;

  // Create arrays to store the positional and
  // orientational data
  darray xc_list(Numcycle*step + 1, 2);
  xc_list = darray::Zero(Numcycle*step + 1, 2);
  darray xc_shift_list(Numcycle*step + 1, 2);
  xc_shift_list = darray::Zero(Numcycle*step + 1, 2);
  dvect theta_list(Numcycle*step + 1);
  theta_list = dvect::Zero(Numcycle*step + 1);

  // Store initial coordinates and theta
  theta_list(0) = theta;
  xc_list(0,0) = xm(0);
  xc_list(0,1) = xm(1);

  dvect x(nvariables);

  std::cout << "Network Construction in Progress..." << std::endl;
  BuildConnections(Gcon, Kcon, BoundaryX, BoundaryY);
  std::cout << "Diluting the network..." << std::endl;
  NetworkDilute(Gcon, 6*pbond);

  std::cout << "Initializing the system..." << std::endl;
  // Initialize all the variables
  CrossLinks(x, r, xm, 0);
  // Shear the network with input (0.0) and
  // compress with eps
  AffineShear(x, 0.0);

  // temporary variables
  double Emin;
  int it;

  // Place the cell on the lattice
  totalnum = adhesion_site(x, xm, theta);

  // initialize the system and
  // consider cell contraction caused by non-zero phase
  contraction(Adhesion, theta, phase);

  return 0;
}
