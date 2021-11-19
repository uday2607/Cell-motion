#ifndef GLOBAL_DEFS
#define GLOBAL_DEFS

#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <random>

// Boost libraries
#include <boost/geometry/core/cs.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/ring.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/segment.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/multi/geometries/multi_point.hpp>
#include <boost/geometry/multi/geometries/multi_linestring.hpp>
#include <boost/geometry/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/geometry/algorithms/intersection.hpp>

// Eigen library
#include <eigen3/Eigen/Core>

//Random Number Generators
#include "../Random/pcg_random.hpp"

// Seed with a real random value, if available
pcg_extras::seed_seq_from<std::random_device> seed_source;

// Make a random number engine
pcg32 rng(seed_source);

// Type definitions
// Eigen typedefs
using Eigen::ArrayXd                                    = dvect;
using Eigen::ArrayXXd                                   = darray;
using Eigen::ArrayXf                                    = fvect;
using Eigen::ArrayXXf                                   = farray;
using Eigen::ArrayXi                                    = ivect;
using Eigen::ArrayXXi                                   = iarray;
using Eigen::Array<bool,Eigen::Dynamic,1>               = bvect;
using Eigen::Array<bool,Eigen::Dynamic,Eigen::Dynamic>  = barray;

// Boost typedefs
using value_type       = double;
using cs_type          = bg::cs::cartesian;
using point_type       = bg::model::point<value_type, 2, cs_type>;
using polygon_type     = bg::model::ring<point_type>;
using line_string_type = bg::model::linestring<point_type>;
using multi_line_type  = bg::model::multi_linestring<line_string_type>;
using intersection     = bg:intersection;

#endif
