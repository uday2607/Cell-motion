#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <vector>
#include <random>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include "../src/random/pcg_random.hpp"


// Seed with a real random value, if available
pcg_extras::seed_seq_from<std::random_device> seed_source;

// Make a random number engine
pcg32 rng(seed_source);

// Typedef
typedef Eigen::ArrayXd dvect;
typedef Eigen::ArrayXXd darray;
typedef Eigen::ArrayXf fvect;
typedef Eigen::ArrayXXf farray;
typedef Eigen::ArrayXi ivect;
typedef Eigen::ArrayXXi iarray;
typedef Eigen::Array<bool,Eigen::Dynamic,1> bvect;
typedef Eigen::Array<bool,Eigen::Dynamic,Eigen::Dynamic> barray;

int main(int argc, char const *argv[]) {

  darray r(3, 2);
  r = darray::Zero(3,2);
  r(1, 1) = 2;
  std::cout << r(2, 1) << std::endl;
  return 0;
}
