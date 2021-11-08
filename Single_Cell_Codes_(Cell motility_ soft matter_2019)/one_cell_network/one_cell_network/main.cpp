#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <string>
#include <fstream>
#include <sstream>
#include "vector.h"
#include "nrutil.h"
#include "dlinmin.h"
#include "dF1dim.h"
#include "frprmn.h"
#include "dbrent.h"
//#include <vector>
#include <random>

//-----global parameters---------------------
int step;
double a;   //long semiaxis
double b;   //short semiaxis
int N;      //number of adhesions
double kp;
double k_b; //back end dissociation rate
double k_f; //front end dissociation rate
double ks;
double dt;  //based on \tau=100dt=1min
double tau;
double Delta;
double lambda;  //strength of contraction
double Alpha;   //dimensionless parameter
int **Gcon;
double **Kcon;

bool **BoundaryX, **BoundaryY;
int Lx, Ly;
double h;
double r[3][2];
int nvariables;
std::string out_dir;
int totalnum;     //total number of adhesion site
//int totalnum2;
double pbond;
double kappa;
double mu;      //stretching stiffness
int Numcycle;   //number of cell cycle to simulate
double fThreshold;
double onrate;
double eps;     //axial extension strain...
int phase1;     //the phase of contraction of cell1
//int phase2;     //the phase of contraction of cell2
//-------------------------------------------

//--------global variables-------------------
double **Adhesion;
//double **Adhesion2;
double **Adhesion0;
//double **Adhesion02;
double **Adhesion_shift;
//double **Adhesion_shift2;
int **AdhesionIndex;
//int **AdhesionIndex2;
int *ifremove;
//int *ifremove2;
//-------------------------------------------

using namespace std;

#include "ifcell.h"
#include "in_out.h"
#include "contraction.h"
#include "extension.h"
#include "energy_grad.h"
#include "rotation.h"
#include "Network.h"
#include "Deformations.h"

double *zerovector(int nl) {
    double *newvector;

    newvector = new double[nl + 1];

    for (int i = 1; i <= nl; i++) {
        newvector[i] = 0.0;
    }

    return newvector;
}


int main(int argc, char *argv[]) {

    out_dir = argv[1];
/*
    ofstream test_file;
    test_file.exceptions(std::ifstream::failbit | std::ifstream::badbit);
    test_file.open(string(argv[1]) + "dyy_test_file.txt");
    test_file << "This is only a test!" << endl;
    test_file.close();
*/

    srand(10000 * atoi(argv[3]));
    cout << "Random seed=" << atoi(argv[3]) << endl;

    //------------------------------------
    step = 30;
    a = 6;
    b = 4;
    N = 4 * (floor(a) + 1) * (floor(b) + 1);
    kp = 1e-2;
    k_b = 2e-2;
    k_f = k_b / 4;
    ks = 0.1;   //original value is 1e-4
    dt = 60.0 / double(step);
    eps = atof(argv[2]) * 0.20;   //normal extension
    tau = 60.0;
    Delta = 0.5e-9;
    lambda = 0.4;
    Alpha = 25;
    fThreshold = 0.005 * ks;
    onrate = 0.4;
    phase1 = 0;       //the first cell always starts with 0 phase;
//    phase2 = rand() % step;     //the second cell starts with random phase;

    Lx = Ly = 60;    //network size
    h = Ly * 0.86602540378443864676 * (1 + eps);   //height after axial extension
    pbond = 1e16;   //bond existence probability
    kappa = 1e16;  //bending stiffness
    mu = 1e16;
    Numcycle = 10;
    totalnum = 0;    //total number of adhesion site
//    totalnum2 = 0;
    //------------------------------------

    nvariables = 1 + 2 * Lx * Ly + 3;

    OutputLx(argv[2], argv[3]);

    //----allocate memory for global variables----
    Adhesion = new double *[N + 1];    //only part of adhesion is utilized
    for (int i = 0; i <= N; i++)
        Adhesion[i] = new double[2];

    Adhesion0 = new double *[N + 1];
    for (int i = 0; i <= N; i++)
        Adhesion0[i] = new double[2];

    Adhesion_shift = new double *[N + 1];
    for (int i = 0; i <= N; i++)
        Adhesion_shift[i] = new double[2];

    AdhesionIndex = new int *[N + 1];
    for (int i = 0; i <= N; i++)
        AdhesionIndex[i] = new int[2];

    ifremove = new int[N + 1];

    //-------------------------second cell--------
/*
    Adhesion2 = new double *[N + 1];    //only part of adhesion is utilized
    for (int i = 0; i <= N; i++)
        Adhesion2[i] = new double[2];

    Adhesion02 = new double *[N + 1];
    for (int i = 0; i <= N; i++)
        Adhesion02[i] = new double[2];

    Adhesion_shift2 = new double *[N + 1];
    for (int i = 0; i <= N; i++)
        Adhesion_shift2[i] = new double[2];

    AdhesionIndex2 = new int *[N + 1];
    for (int i = 0; i <= N; i++)
        AdhesionIndex2[i] = new int[2];

    ifremove2 = new int[N + 1];
*/

    BoundaryX = new bool *[3];
    for (int i = 0; i <= 2; i++)
        BoundaryX[i] = new bool[1 + Lx * Ly];

    BoundaryY = new bool *[3];
    for (int i = 0; i <= 2; i++)
        BoundaryY[i] = new bool[1 + Lx * Ly];

    Gcon = new int *[3];
    for (int i = 0; i <= 2; i++)
        Gcon[i] = new int[1 + Lx * Ly];

    Kcon = new double *[3];
    for (int i = 0; i <= 2; i++)
        Kcon[i] = new double[1 + Lx * Ly];
    //--------------------------------------------

    //----displacement vector------------xnew=xold+r[con][0], ynew=yold+r[con][1]
    r[0][0] = -0.5;
    r[0][1] = 0.86602540378443864676;

    r[1][0] = 0.5;
    r[1][1] = 0.86602540378443864676;

    r[2][0] = 1;
    r[2][1] = 0;
    //----------------------------------------------

    //-----declare local variables----------------
    double xm[2];  //center of the cell
    xm[0] = double(Lx) * (3.0 / 4.0) - 10;
    xm[1] = double(Ly) / 2 * sqrt(3) / 2;

/*
    double xm2[2];
    xm2[0] = double(Lx) * (3.0 / 4.0) + 10;
    xm2[1] = double(Ly) / 2 * sqrt(3) / 2;
*/

//    double theta=distribution(generator);
    double theta=fRand(0., 3.14159265358979 * 2.);
//    double theta = 1.57;
//    double theta=0.0;

    cout << "theta0= " << theta << endl;

    double theta0;            //the orientation of the cell
    double xc_list[Numcycle * step + 1][2];
    double theta_list[Numcycle * step + 1];
    double xc_shift_list[Numcycle * step + 1][2];

/*
    double theta2=fRand(0, 6.28);
    double theta2 = 1.57;
    double theta2=0.0;
    cout << "theta02= " << theta2 << endl;

    double theta02;
    double xc_list2[Numcycle * step + 1][2];
    double theta_list2[Numcycle * step + 1];
*/

    for (int i = 0; i <= Numcycle * step; i++) {
        xc_list[i][0] = 0.0;
        xc_list[i][1] = 0.0;
        xc_shift_list[i][0] = 0.;
        xc_shift_list[i][1] = 0.;
        theta_list[i] = 0.0;
/*
        xc_list2[i][0] = 0.0;
        xc_list2[i][1] = 0.0;
        theta_list2[i] = 0.0;
*/
    }
    theta_list[0] = theta;
    xc_list[0][0] = xm[0];
    xc_list[0][1] = xm[1];

/*
    xc_list2[0][0] = xm2[0];
    xc_list2[0][1] = xm2[1];
*/

    double *x;
    x = new double[nvariables];

    cout << "network construction" << endl;
    ConnectionsBuildFull(Gcon);
    cout << "dilute the network" << endl;
    Dilute(Gcon, 6 * pbond);


    cout << "initialize network" << endl;
    CrossLinks(x, xm, 0);   //initialize all variables  dtheta1=0, dtheta2=0;
    AffineShear(x, 0.0);   //shear with input (0.0) and compress with eps the network

    double Emin;
    int it;

    //put in the cells
    totalnum = adhesion_site(x, xm, theta);
//    totalnum2 = adhesion_site2(x, xm2, theta2);

//    mature(x, theta, xm[0], xm[1]);
//    mature2(x, theta2, xm2[0], xm2[1]);

    //initialize the system, consider cell contraction caused by non-zero phase
    contraction(Adhesion, theta, phase1);
//    contraction2(Adhesion2, theta2, phase2);

    frprmn(x, nvariables - 1, 1e-8, &it, &Emin, energy, grad);

//    for (int nstep = 1; nstep <= 3; nstep++) {
    for (int nstep = 1; nstep <= Numcycle * step; nstep++) {

        cout << "----------------------------------" << endl;
        cout << "nstep=" << nstep << endl;
        cout << "----------------------------------" << endl;

        contraction(Adhesion, theta, 1);     //contract the cell
//        contraction2(Adhesion2, theta2, 1);

        int iter;
        double EnergyMin;

        x[nvariables - 1] = 0.0;    //initialize dtheta1 and dtheta2 for better optimization
//        x[nvariables - 4] = 0.0;

        cout << "step=" << nstep << endl;
        frprmn(x, nvariables - 1, 1e-8, &iter, &EnergyMin, energy, grad);

        cout << "-------------------------" << endl;
        cout << x[nvariables - 3] - xm[0] << "," << x[nvariables - 2] - xm[1] << "," << x[nvariables - 1] << endl;
//        cout << x[nvariables - 6] - xm[0] << "," << x[nvariables - 5] - xm[1] << "," << x[nvariables - 4] << endl;

        //---save the results-------------------
/*
        xc_list[nstep][0] = x[nvariables - 6];           //x
        xc_list[nstep][1] = x[nvariables - 5];           //y

        xc_list2[nstep][0] = x[nvariables - 3];           //x
        xc_list2[nstep][1] = x[nvariables - 2];           //y
*/

        xc_list[nstep][0] = x[nvariables - 3];
        xc_list[nstep][1] = x[nvariables - 2];

//        rotation(x[nvariables - 4]);      //rotate the cell
//        rotation2(x[nvariables - 1]);
        rotation(x[nvariables - 1]);
//stops here
        // update the angle
/*
        theta = theta + x[nvariables - 4];
        theta_list[nstep] = theta;

        theta2 = theta2 + x[nvariables - 1];
        theta_list2[nstep] = theta2;
*/

        theta = theta + x[nvariables - 1];
        theta_list[nstep] = theta;

//save new net config & cell location & orientation---
//        OutputAd(Adhesion, ifremove, nstep, totalnum);
//
//        OutputAd2(Adhesion2, ifremove2, nstep, totalnum2);
//
//        Outputnodes(x, nstep, argv[2], argv[3]);

        phase1 = (phase1 + 1) % step;
//        phase2 = (phase2 + 1) % step;

        if (phase1 == 0) {

            double xc, yc;
            xc = x[nvariables - 3];
            yc = x[nvariables - 2];
            extension(Adhesion, xc, yc, xm, theta);  //shift the cell center, relax the cell
            double xm_before_shift[2] = {xm[0], xm[1]};
            shiftcell(x, xm, theta);
            xc_shift_list[nstep][0] = xm[0] - xm_before_shift[0];
            xc_shift_list[nstep][1] = xm[1] - xm_before_shift[1];
            theta0 = theta;
            totalnum = adhesion_site(x, xm, theta);
        } else if (phase1 == 1) {
            mature(x, theta0, xm[0], xm[1]);
        } else {
            detach(x, theta0, xm[0], xm[1]);
        }

/*
        if (phase2 == 0) {

            double xc2, yc2;
            xc2 = x[nvariables - 3];
            yc2 = x[nvariables - 2];
            extension2(Adhesion2, xc2, yc2, xm2, theta2);

            shiftcell2(x, xm2, theta2);
            theta02 = theta2;
            totalnum2 = adhesion_site2(x, xm2, theta2);
        } else if (phase2 == 1) {
            mature2(x, theta02, xm2[0], xm2[1]);
        } else {
            detach2(x, theta02, xm2[0], xm2[1]);
        }
*/


    }

//    cout << "phase2=" << phase2 << endl;

    OutputXc(xc_list, theta_list, argv[2], argv[3]);
//    OutputXc2(xc_list2, theta_list2, argv[2], argv[3]);

    OutputShift(xc_shift_list, argv[2], argv[3]);

    OutputPhase(phase1, "phase.txt", argv[2], argv[3]);
//    Outputbonds(x, "bonds.txt",argv[2], argv[3]);


//    cout<<"stretching="<<EnergyStretch(x)<<endl;
//    cout<<"bending="<<kappa*EnergyBend(x)<<endl;
//    cout<<"ad1="<<AdEnergy(x)<<endl;
//    cout<<"ad2="<<AdEnergy2(x)<<endl;

    return 0;
}
