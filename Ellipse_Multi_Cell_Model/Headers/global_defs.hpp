#ifndef GLOBAL_DEFS
#define GLOBAL_DEFS

#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <random>
#include <eigen3/Eigen/Core>

//Random Number Generators
#include "../Random/pcg_random.hpp"

// Seed with a real random value, if available
pcg_extras::seed_seq_from<std::random_device> seed_source;

// Make a random number engine
pcg32 rng(seed_source);

// Type definitions
typedef Eigen::ArrayXd dvect;
typedef Eigen::ArrayXXd darray;
typedef Eigen::ArrayXf fvect;
typedef Eigen::ArrayXXf farray;
typedef Eigen::ArrayXi ivect;
typedef Eigen::ArrayXXi iarray;
typedef Eigen::Array<bool,Eigen::Dynamic,1> bvect;
typedef Eigen::Array<bool,Eigen::Dynamic,Eigen::Dynamic> barray;

#endif
