#ifndef UTILS
#define UTILS

#include "global_defs.hpp"

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

// Comparision functions
class isGreaterDouble{
  public:
    vector<int>* GT;
    double G;
    isGreater(vector<int> *g, double a){GT = g; G = a;}
    void operator()(double i){static int it = 0;
                    if (i > G) GT->push_back(it); it++;}
};

class isLesserDouble{
  public:
    vector<int>* LT;
    double L;
    isLesser(vector<int> *l, double a){LT = l; L = a;}
    void operator()(double i){static int it = 0;
                    if (i < L) LT->push_back(it); it++;}
};

class isEqualDouble{
  public:
    vector<int>* ET;
    double E;
    isEqual(vector<int> *e, double a){ET = e; E = a;}
    void operator()(double i){static int it = 0;
                    if (i == E) ET->push_back(it); it++;}
};

class isGreaterInt{
  public:
    vector<int>* GT;
    int G;
    isGreater(vector<int> *g, int a){GT = g; G = a;}
    void operator()(int i){static int it = 0;
                    if (i > G) GT->push_back(it); it++;}
};

class isLesserInt{
  public:
    vector<int>* LT;
    int L;
    isLesser(vector<int> *l, int a){LT = l; L = a;}
    void operator()(int i){static int it = 0;
                    if (i < L) LT->push_back(it); it++;}
};

class isEqualInt{
  public:
    vector<int>* ET;
    int E;
    isEqual(vector<int> *e, int a){ET = e; E = a;}
    void operator()(int i){static int it = 0;
                    if (i == E) ET->push_back(it); it++;}
};

#endif
