// Declaration of global variables
#ifndef GLOBAL_H
#define GLOBAL_H

#include <string>
#include <eigen3/Eigen/Core>

typedef Eigen::ArrayXd dvect;
typedef Eigen::ArrayXXd darray;
typedef Eigen::ArrayXf fvect;
typedef Eigen::ArrayXXf farray;
typedef Eigen::ArrayXi ivect;
typedef Eigen::ArrayXXi iarray;
typedef Eigen::Array<bool,Eigen::Dynamic,1> bvect;
typedef Eigen::Array<bool,Eigen::Dynamic,Eigen::Dynamic> barray;

// Global Parameters
//-time params-
extern int step;
extern double dt;           //based on tau = 100dt = 1min
extern double tau;          //contraction cycle time scale

//-cell params-
extern double a;            //long semiaxis
extern double b;            //short semiaxis
extern int N;               //number of adhesions
extern double kp;           //onrate parameter
extern double onrate;       //onrate
extern double k_b;          //back end dissociation rate
extern double k_f;          //front end dissociation rate
extern double ks;           //FA spring constant
extern double lambda;       //strength of contraction
extern double Alpha;        //strength of the bond
extern double Delta;        //
extern int totalnum;        //total number of adhesion sites
extern int phase;           //the phase of contraction of cell
                            //initial = 0
//-lattice-params-
extern double eps;          //axial extension strain
                            //user input
extern double pbond;        //bond existence probability
extern double kappa;        //bending stiffness
extern double mu;           //stretching stiffness
extern double fThreshold;   //maturation force threshold
extern int Lx;
extern int Ly;              //network size
extern double h;            //height after axial extension

//-numerical-stuff-
extern int Numcycle;        //number of cell cycle to simulate
extern int nvariables;      //total number of variables
extern std::string out_dir;

// Global Variables
extern darray Adhesion;
extern darray Adhesion0;
extern darray Adhesion_shift;
extern iarray AdhesionIndex;
extern ivect ifremove;
extern iarray Gcon;
extern darray Kcon;
extern barray BoundaryX;
extern barray BoundaryY;
extern darray r;

// Definitions
// Global Parameters
//-time params-
int step = 30;
double dt = 60/double(step);           //based on tau = 100dt = 1min
double tau = 60.0;                     //contraction cycle time scale

//-cell params-
double a = 6;                          //long semiaxis
double b = 4;                          //short semiaxis
int N = 4*(floor(a)+1)*(floor(b)+1);   //number of adhesions
double kp = 1e-2;                      //onrate parameter
double onrate = 0.4;                   //onrate
double k_b = 2e-2;                     //back end dissociation rate
double k_f = 5e-3;                     //front end dissociation rate
double ks = 0.1;                       //FA spring constant
double lambda = 0.4;                   //strength of contraction
double Alpha = 25;                     //strength of the bond
double Delta = 0.5e-9;                 //
int totalnum = 0;                      //total number of adhesion sites
int phase = 0;                         //the phase of contraction of cell
                                       //initial = 0
//-lattice-params-
double eps;                            //axial extension strain
                                       //user input
double pbond = 1.0;                    //bond existence probability
double kappa = 1e16;                   //bending stiffness
double mu = 1;                         //stretching stiffness
double fThreshold = 0.005*ks;          //maturation force threshold
int Lx = 60;
int Ly = 60;                           //network size
double h;                              //height after axial extension

//-numerical-stuff-
int Numcycle = 10;                     //number of cell cycle to simulate
int nvariables = 1 + 2*Lx*Ly + 3;      //total number of variables
std::string out_dir;

// Global Variables
darray Adhesion(N+1, 2);
darray Adhesion0(N+1, 2);
darray Adhesion_shift(N+1, 2);
iarray AdhesionIndex(N+1, 2);
ivect ifremove(N+1);
iarray Gcon(3, 1+Lx*Ly);
darray Kcon(3, 1+Lx*Ly);
barray BoundaryX(3, 1+Lx*Ly);
barray BoundaryY(3, 1+Lx*Ly);
darray r(3, 2);

#endif
