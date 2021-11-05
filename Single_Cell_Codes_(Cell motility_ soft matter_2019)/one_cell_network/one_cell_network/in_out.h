void OutputAd(double **Adhesion, int *ifremove, int step, int totalnum) {
    ofstream Ainfo;
    Ainfo.open(out_dir + "Adhesion_step" + std::to_string(step) + ".txt");

    for (int num = 1; num <= totalnum; num++) {
        Ainfo << Adhesion[num][0] << " " << Adhesion[num][1] << " " << AdhesionIndex[num][0] << " "
              << AdhesionIndex[num][1] << " " << Adhesion_shift[num][0] << " " << Adhesion_shift[num][1] << " "
              << ifremove[num] << endl;
    }

    Ainfo.close();
}

/*void OutputAd2(double **Adhesion2, int *ifremove2, int step, int totalnum2) {
    ofstream Ainfo;
    Ainfo.open(out_dir + "Adhesion2_step" + std::to_string(step) + ".txt");

    for (int num = 1; num <= totalnum2; num++) {
        Ainfo << Adhesion2[num][0] << " " << Adhesion2[num][1] << " " << AdhesionIndex2[num][0] << " "
              << AdhesionIndex2[num][1] << " " << Adhesion_shift2[num][0] << " " << Adhesion_shift2[num][1] << " "
              << ifremove2[num] << endl;
    }

    Ainfo.close();
}*/

void OutputXc(double xc_list[][2], double theta_list[], char *frac, char *index) {
    ofstream xcinfo;
    ofstream thetainfo;

    xcinfo.open(out_dir + "f" + frac + "id" + index + "xc.txt");
    thetainfo.open(out_dir + "f" + frac + "id" + index + "theta.txt");

    for (int num = 0; num <= Numcycle * step; num++) {
        xcinfo << xc_list[num][0] << " " << xc_list[num][1] << endl;
        thetainfo << theta_list[num] << endl;
    }

    xcinfo.close();
    thetainfo.close();
}

void OutputShift(double shift_list[][2], char *frac, char *index) {
    ofstream shift_info;

    shift_info.open(out_dir + "f" + frac + "id" + index + "xcShift.txt");

    for (int num = 0; num <= Numcycle * step; num++) {
        shift_info << shift_list[num][0] << " " << shift_list[num][1] << endl;
    }

    shift_info.close();
}

/*
void OutputXc2(double xc_list2[901][2], double theta_list2[901], char *frac, char *index) {
    ofstream xcinfo;
    ofstream thetainfo;

    xcinfo.open(out_dir + "xc2" + "f" + frac + "id" + index + ".txt");
    thetainfo.open(out_dir + "theta2" + "f" + frac + "id" + index + ".txt");

    for (int num = 1; num <= Numcycle * step; num++) {
        xcinfo << xc_list2[num][0] << " " << xc_list2[num][1] << endl;
        thetainfo << theta_list2[num] << endl;
    }

    xcinfo.close();
    thetainfo.close();
}
*/

void Outputnodes(double *x, int counter, char *frac, char *index) {
    ofstream nodeinfo;
//    cout<<out_dir+"nodes"+std::to_string(counter)+"f"+frac+"id"+index+".txt"<<endl;
    nodeinfo.open(out_dir + std::to_string(counter) + "f" + frac + "id" + index + "nodes.txt");

    int N = Lx * Ly, G1;
    for (int i = 0; i <= (Lx - 1); i++) {
        for (int j = 0; j <= (Ly - 1); j++) {
            G1 = 1 + j * Lx + i;
            nodeinfo << x[G1] << " " << x[G1 + N] << " " << G1 << endl;
        }
    }

    nodeinfo.close();
}

void OutputLx(string frac, string job_index) {
    ofstream info;
    info.exceptions(std::ofstream::failbit | std::ofstream::badbit);

    info.open(out_dir + "f" + frac + "id" + job_index + "Lx.txt");
    info << Lx << endl;

    info.close();
}

void OutputPhase(int phase, string s1, char *frac, char *index) {
    ofstream phaseinfo;
    phaseinfo.open(out_dir + "f" + frac + "id" + index + s1);

    phaseinfo << phase << endl;

    phaseinfo.close();
}

void Outputbonds(double *x, string s1, char *frac, char *index) {
    ofstream bondinfo;
    bondinfo.open(out_dir + "f" + frac + "id" + index + s1);

    int N = Lx * Ly, G1, G2;
    for (int i = 0; i <= (Lx - 1); i++) {
        for (int j = 0; j <= (Ly - 1); j++) {
            G1 = 1 + j * Lx + i;

            for (int con = 0; con <= 2; con++) {
                G2 = Gcon[con][G1];

                if (G2 != -1 && !BoundaryX[con][G1] && !BoundaryY[con][G1]) {

                    bondinfo << G1 << " " << G2 << " " << endl;
                }

            }

        }
    }

    bondinfo.close();
}
