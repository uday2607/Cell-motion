double fRand(double fMin, double fMax) {
    if (fMin == fMax) {
        return fMin;
    } else {
        double f = (double) rand() / (double) RAND_MAX;
        return fMin + f * (fMax - fMin);
    }
}

void rotation(double dtheta) {

    for (int i = 1; i <= N; i++) {  //should be totalnum, but the results not affected
        double x = Adhesion[i][0];
        double y = Adhesion[i][1];

        Adhesion[i][0] = cos(dtheta) * x - sin(dtheta) * y;
        Adhesion[i][1] = sin(dtheta) * x + cos(dtheta) * y;
    }

}

/*void rotation2(double dtheta) {

    for (int i = 1; i <= N; i++) {
        double x = Adhesion2[i][0];
        double y = Adhesion2[i][1];

        Adhesion2[i][0] = cos(dtheta) * x - sin(dtheta) * y;
        Adhesion2[i][1] = sin(dtheta) * x + cos(dtheta) * y;
    }

}*/

void mature(double *x, double theta0, double xc0, double yc0) {

//    double xc = x[nvariables - 6];
//    double yc = x[nvariables - 5];

    double xc = x[nvariables - 3];
    double yc = x[nvariables - 2];

    cout << "--------------------------------------" << endl;

    for (int i = 1; i <= totalnum; i++) {

        if (ifremove[i] == 1) continue;

        double x0 = (Adhesion0[i][0] - xc0) * cos(theta0) +
                    (Adhesion0[i][1] - yc0) * sin(theta0);  //make sure this is correct
        double y0 = -(Adhesion0[i][0] - xc0) * sin(theta0) + (Adhesion0[i][1] - yc0) * cos(theta0);
        double x1 = x[AdhesionIndex[i][0]] + Adhesion_shift[i][0];
        double y1 = x[AdhesionIndex[i][1]] + Adhesion_shift[i][1];

        double k_x;

        k_x = k_b - (k_b - k_f) * (x0 + a) / (2 * a);   //basal detaching rate is determined by the initial position

        double xx = Adhesion[i][0] + xc;
        double yy = Adhesion[i][1] + yc;

        double dis = sqrt((xx - x1) * (xx - x1) + (yy - y1) * (yy - y1));
        double force = ks * dis;
        double off_rate = exp(-force / fThreshold);

//        cout<<"i="<<i<<",rx= "<<(Adhesion0[i][0]-xc0)<<", ry= "<< (Adhesion0[i][1]-yc0)<<", rx2= "<<x0<<", ry2="<<y0<<endl;

        if (fRand(0, 1) < off_rate * dt) {
            ifremove[i] = 1;        //the adhesion site is detached
            cout << "one nascent adhesion site perishes" << endl;
        }
    }
}

/*void mature2(double *x, double theta02, double xc02, double yc02) {

    double xc = x[nvariables - 3];
    double yc = x[nvariables - 2];

    cout << "--------------------------------------" << endl;

    for (int i = 1; i <= totalnum2; i++) {

        if (ifremove2[i] == 1) continue;

        double x0 = (Adhesion02[i][0] - xc02) * cos(theta02) +
                    (Adhesion02[i][1] - yc02) * sin(theta02);  //make sure this is correct
        double y0 = -(Adhesion02[i][0] - xc02) * sin(theta02) + (Adhesion02[i][1] - yc02) * cos(theta02);
        double x1 = x[AdhesionIndex2[i][0]] + Adhesion_shift2[i][0];
        double y1 = x[AdhesionIndex2[i][1]] + Adhesion_shift2[i][1];

        double k_x;

        k_x = k_b - (k_b - k_f) * (x0 + a) / (2 * a);   //basal detaching rate is determined by the initial position

        double xx = Adhesion2[i][0] + xc;
        double yy = Adhesion2[i][1] + yc;

        double dis = sqrt((xx - x1) * (xx - x1) + (yy - y1) * (yy - y1));
        double force = ks * dis;
        double off_rate = exp(-force / fThreshold);

//        cout<<"i="<<i<<",rx= "<<(Adhesion0[i][0]-xc0)<<", ry= "<< (Adhesion0[i][1]-yc0)<<", rx2= "<<x0<<", ry2="<<y0<<endl;

        if (fRand(0, 1) < off_rate * dt) {
            ifremove2[i] = 1;        //the adhesion site is detached
            cout << "one adhesion site breaks" << endl;
        }
    }
}
*/


bool ismature() {

    int count = 0;

    for (int i = 1; i <= totalnum; i++) {
        if (ifremove[i] == 1) count++;
    }

    return count != totalnum;
}

/*
bool ismature2() {

    int count = 0;

    for (int i = 1; i <= totalnum2; i++) {
        if (ifremove2[i] == 1) count++;
    }

    if (count == totalnum2) return false;
    else return true;
}
*/


void detach(double *x, double theta0, double xc0, double yc0) {


//    cout<<"xc0= "<<xc0<<",yc0="<<yc0<<", theta0="<<theta0<<endl;

//    double xc = x[nvariables - 6];
//    double yc = x[nvariables - 5];

    double xc = x[nvariables - 3];
    double yc = x[nvariables - 2];

    cout << "--------------------------------------" << endl;

    for (int i = 1; i <= totalnum; i++) {

        if (ifremove[i] == 1) continue;

        double x0 = (Adhesion0[i][0] - xc0) * cos(theta0) +
                    (Adhesion0[i][1] - yc0) * sin(theta0);  //make sure this is correct
        double y0 = -(Adhesion0[i][0] - xc0) * sin(theta0) + (Adhesion0[i][1] - yc0) * cos(theta0);
        double x1 = x[AdhesionIndex[i][0]] + Adhesion_shift[i][0];
        double y1 = x[AdhesionIndex[i][1]] + Adhesion_shift[i][1];

        double k_x;

        k_x = k_b - (k_b - k_f) * (x0 + a) / (2 * a);   //basal detaching rate is determined by the initial position

        double xx = Adhesion[i][0] + xc;
        double yy = Adhesion[i][1] + yc;

        double dis = sqrt((xx - x1) * (xx - x1) + (yy - y1) * (yy - y1));
        double off_rate = k_x * exp(Alpha * dis / a * ks / 0.1);

//        cout<<"i="<<i<<",rx= "<<(Adhesion0[i][0]-xc0)<<", ry= "<< (Adhesion0[i][1]-yc0)<<", rx2= "<<x0<<", ry2="<<y0<<endl;

        if (fRand(0, 1) < off_rate * dt) {
            ifremove[i] = 1;        //the adhesion site is detached
            cout << "one adhesion site breaks" << endl;
        }
    }

}

/*
void detach2(double *x, double theta02, double xc02, double yc02) {


//    cout<<"xc0= "<<xc0<<",yc0="<<yc0<<", theta0="<<theta0<<endl;

    double xc = x[nvariables - 3];
    double yc = x[nvariables - 2];

    cout << "--------------------------------------" << endl;

    for (int i = 1; i <= totalnum2; i++) {

        if (ifremove2[i] == 1) continue;

        double x0 = (Adhesion02[i][0] - xc02) * cos(theta02) +
                    (Adhesion02[i][1] - yc02) * sin(theta02);  //make sure this is correct
        double y0 = -(Adhesion02[i][0] - xc02) * sin(theta02) + (Adhesion02[i][1] - yc02) * cos(theta02);
        double x1 = x[AdhesionIndex2[i][0]] + Adhesion_shift2[i][0];
        double y1 = x[AdhesionIndex2[i][1]] + Adhesion_shift2[i][1];

        double k_x;

        k_x = k_b - (k_b - k_f) * (x0 + a) / (2 * a);   //basal detaching rate is determined by the initial position

        double xx = Adhesion2[i][0] + xc;
        double yy = Adhesion2[i][1] + yc;

        double dis = sqrt((xx - x1) * (xx - x1) + (yy - y1) * (yy - y1));
        double off_rate = k_x * exp(Alpha * dis / a * ks / 0.1);

//        cout<<"i="<<i<<",rx= "<<(Adhesion0[i][0]-xc0)<<", ry= "<< (Adhesion0[i][1]-yc0)<<", rx2= "<<x0<<", ry2="<<y0<<endl;

        if (fRand(0, 1) < off_rate * dt) {
            ifremove2[i] = 1;        //the adhesion site is detached
            cout << "one adhesion site breaks" << endl;
        }
    }

}
*/
