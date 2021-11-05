void contraction(double **Adhesion, double theta, int fold) {

    for (int i = 1; i <= N; i++) {
        double x;
        double y;
        double c0 = lambda * dt / tau * ((double) fold);
//        cout<<"c0="<<c0<<endl;

        x = Adhesion[i][0];
        y = Adhesion[i][1];

        Adhesion[i][0] = x - c0 * cos(theta) * cos(theta) * x - c0 * sin(theta) * cos(theta) * y;
        Adhesion[i][1] = y - c0 * sin(theta) * cos(theta) * x - c0 * sin(theta) * sin(theta) * y;

    }
}

/*
void contraction2(double **Adhesion2, double theta2, int fold) {

    for (int i = 1; i <= N; i++) {
        double x;
        double y;
        double c0 = lambda * dt / tau * ((double) fold);
//        cout<<"c0="<<c0<<endl;

        x = Adhesion2[i][0];
        y = Adhesion2[i][1];

        Adhesion2[i][0] = x - c0 * cos(theta2) * cos(theta2) * x - c0 * sin(theta2) * cos(theta2) * y;
        Adhesion2[i][1] = y - c0 * sin(theta2) * cos(theta2) * x - c0 * sin(theta2) * sin(theta2) * y;

    }
}
*/
