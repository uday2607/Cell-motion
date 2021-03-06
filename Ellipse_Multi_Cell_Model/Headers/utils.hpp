#ifndef UTILS
#define UTILS

#include "global_defs.hpp"
using namespace std;
using namespace Eigen;

// Global Random functions
// Uniform distribution
double randDouble(double low, double high) {
  std::uniform_real_distribution<double> dist(low, high);
  return dist(rng);
}

// Uniform integer distribution
int randInt(int low, int high) {
  std::uniform_int_distribution<int> dist(low, high);
  return dist(rng);
}

// Round off Function
double round_up(double value, int decimal_places) {
    const double multiplier = std::pow(10.0, decimal_places);
    return std::ceil(value * multiplier) / multiplier;
}

// Search functions for Eigen vectors
int Search_dvect_int(dvect vect, int x) {

  for (int i = 0; i < vect.size(); i++) {
    if (vect(i) == x) {
      return i;
    }
  }

  return -1;
}

int Search_ivect_int(ivect vect, int x) {

  for (int i = 0; i < vect.size(); i++) {
    if (vect(i) == x) {
      return i;
    }
  }

  return -1;
}

// Comparision functions
class isGreaterDouble{
  public:
    vector<int>* GT;
    double G;
    isGreaterDouble(vector<int> *g, double a){GT = g; G = a;}
    void operator()(double i){static int it = 0;
                    if (i > G) GT->push_back(it); it++;}
};

class isLesserDouble{
  public:
    vector<int>* LT;
    double L;
    isLesserDouble(vector<int> *l, double a){LT = l; L = a;}
    void operator()(double i){static int it = 0;
                    if (i < L) LT->push_back(it); it++;}
};

class isEqualDouble{
  public:
    vector<int>* ET;
    double E;
    isEqualDouble(vector<int> *e, double a){ET = e; E = a;}
    void operator()(double i){static int it = 0;
                    if (i == E) ET->push_back(it); it++;}
};

class isGreaterInt{
  public:
    vector<int>* GT;
    int G;
    isGreaterInt(vector<int> *g, int a){GT = g; G = a;}
    void operator()(int i){static int it = 0;
                    if (i > G) GT->push_back(it); it++;}
};

class isLesserInt{
  public:
    vector<int>* LT;
    int L;
    isLesserInt(vector<int> *l, int a){LT = l; L = a;}
    void operator()(int i){static int it = 0;
                    if (i < L) LT->push_back(it); it++;}
};

class isEqualInt{
  public:
    vector<int>* ET;
    int E;
    isEqualInt(vector<int> *e, int a){ET = e; E = a;}
    void operator()(int i){static int it = 0;
                    if (i == E) ET->push_back(it); it++;}
};

// Numpy functions
void dlinspace(dvect &vect, double low, double high, int num) {
  for (int i = 0; i < num; i++) {
    vect(i) = low + i*(high-low)/N;
  }
}

#endif
