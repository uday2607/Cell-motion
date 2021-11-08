// Standard Libraries
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <fstream>
#include <vector>
#include <random>
#include <eigen3/Eigen/Core>
#include "random/pcg_random.hpp"

// Custom Libraries
#include "Headers/globals.hpp"
#include "Headers/fileio.hpp"
//#include "Headers/network.hpp"


// Seed with a real random value, if available
pcg_extras::seed_seq_from<std::random_device> seed_source;

// Make a random number engine
pcg32 rng(seed_source);

// Uniform distribution
double randDouble(double low, double high) {
  std::uniform_real_distribution<double> dist(low, high);
  return dist(rng);
}

int main(int argc, char const *argv[]) {

  //Output
  out_dir = "Output/" + std::string(argv[1]);

  // Store Length of the lattice used
  std::cout << Ly << std::endl;
  OutputLx(argv[2], argv[3]);

  // Normal extension
  eps = atof(argv[2]) * 0.20;
  h = Ly*0.86602540378443864676*(1+eps);


  // Displacement vector
  // xnew=xold+r[con][0], ynew=yold+r[con][1]
  r(0,0) = -0.5;
  r(0,1) = 0.86602540378443864676;

  r(1,0) = 0.5;
  r(1,1) = 0.86602540378443864676;

  r(2,0) = 1;
  r(2,1) = 0;

  // Local Variables
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

  std::cout << "Network Construction in Progres..." << std::endl;


  return 0;
}
