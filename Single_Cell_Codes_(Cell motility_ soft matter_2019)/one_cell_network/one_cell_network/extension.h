void extension(double **Adhesion, double xc, double yc, double xm[], double theta) {
    //move the cell to the leftmost adhesion site

    double **tempAd;
    tempAd = new double *[totalnum + 1];
    for (int i = 0; i <= totalnum; i++)
        tempAd[i] = new double[2];

    double **tempAdshift;
    tempAdshift = new double *[totalnum + 1];
    for (int i = 0; i <= totalnum; i++)
        tempAdshift[i] = new double[2];

    double *xprime;
    xprime = new double[totalnum + 1];

    int **tempAdIndex;
    tempAdIndex = new int *[totalnum + 1];
    for (int i = 0; i <= totalnum; i++)
        tempAdIndex[i] = new int[2];

    //tempvector is a (totalnum+1) by 2 vector.

    int counter = 1;


    for (int i = 1; i <= totalnum; i++) {

        if (ifremove[i] == 1) continue;  //ignore the site that's already been deleted

        double x;
        double y;

        x = Adhesion[i][0];
        y = Adhesion[i][1];

        tempAd[counter][0] = x;
        tempAd[counter][1] = y;

        //----rotate the reference by theta, should be equilavent to rotate adhesion by -theta------------------
        xprime[counter] = cos(theta) * x + sin(theta) * y;
        //------------------------------------------

        tempAdIndex[counter][0] = AdhesionIndex[i][0];
        tempAdIndex[counter][1] = AdhesionIndex[i][1];

        tempAdshift[counter][0] = Adhesion_shift[i][0];
        tempAdshift[counter][1] = Adhesion_shift[i][1];

        counter++;
    }

    double minval = 1000;
    int minIndex = -1;

    //find the leftmost point
    for (int i = 1; i < counter; i++) {
        if (xprime[i] < minval) {
            minval = xprime[i];
            minIndex = i;
        }
    }

    cout << "minval= " << minval << endl;

    double xm_new[2];

    if (minIndex == -1) {
        xm_new[0] = xc;
        xm_new[1] = yc;
    } else {
        xm_new[0] = xc + 0.8 * (a + minval) * cos(theta);
        xm_new[1] = yc + 0.8 * (a + minval) * sin(theta);
    }

    for (int i = 1; i < counter; i++) {
        tempAd[i][0] = tempAd[i][0] + xc - xm_new[0];
        tempAd[i][1] = tempAd[i][1] + yc - xm_new[1];
    }

    for (int i = 1; i <= totalnum; i++) {

        if (i < counter) {
            Adhesion[i][0] = tempAd[i][0];
            Adhesion[i][1] = tempAd[i][1];
            AdhesionIndex[i][0] = tempAdIndex[i][0];
            AdhesionIndex[i][1] = tempAdIndex[i][1];
            Adhesion_shift[i][0] = tempAdshift[i][0];
            Adhesion_shift[i][1] = tempAdshift[i][1];
        } else {
            Adhesion[i][0] = 0.0;
            Adhesion[i][1] = 0.0;
            AdhesionIndex[i][0] = 0;
            AdhesionIndex[i][1] = 0;
            Adhesion_shift[i][0] = 0;
            Adhesion_shift[i][1] = 0;
        }

        ifremove[i] = 0;
        Adhesion0[i][0] = 0;
        Adhesion0[i][1] = 0;
    }

    xm[0] = xm_new[0];
    xm[1] = xm_new[1];

    //reset the totalnum
    totalnum = counter - 1;

    //free the tempertory pointer
    delete[] xprime;

    for (int i = 0; i <= totalnum; i++)
        delete[] tempAd[i];
    delete[] tempAd;

    for (int i = 0; i <= totalnum; i++)
        delete[] tempAdIndex[i];
    delete[] tempAdIndex;

    for (int i = 0; i <= totalnum; i++)
        delete[] tempAdshift[i];
    delete[] tempAdshift;
    //--------------------------
}


/*
void extension2(double **Adhesion2, double xc2, double yc2, double xm2[], double theta2){
    //move the cell to the leftmost adhesion site

    double **tempAd;
    tempAd = new double*[totalnum2+1];
    for (int i=0; i<=totalnum2; i++)
        tempAd[i]=new double[2];

    double **tempAdshift;
    tempAdshift = new double*[totalnum2+1];
    for (int i=0; i<=totalnum2; i++)
        tempAdshift[i]=new double[2];

    double *xprime;
    xprime = new double [totalnum2+1];

    int **tempAdIndex;
    tempAdIndex= new int*[totalnum2+1];
    for (int i=0; i<=totalnum2; i++)
        tempAdIndex[i]=new int[2];

    //tempvector is a (totalnum+1) by 2 vector.

    int counter=1;


//    cout<<"theta2="<<theta2<<endl;

    for (int i=1; i<=totalnum2; i++){

        if (ifremove2[i]==1) continue;  //ignore the site that's already been deleted

        double x;
        double y;

        x=Adhesion2[i][0];
        y=Adhesion2[i][1];

        tempAd[counter][0]=x;
        tempAd[counter][1]=y;

        //----rotate the reference by theta, should be equilavent to rotate adhesion by -theta------------------
        xprime[counter]=cos(theta2)*x+sin(theta2)*y;
        //------------------------------------------

        tempAdIndex[counter][0]=AdhesionIndex2[i][0];
        tempAdIndex[counter][1]=AdhesionIndex2[i][1];

        tempAdshift[counter][0]=Adhesion_shift2[i][0];
        tempAdshift[counter][1]=Adhesion_shift2[i][1];


//        cout<<"index="<<i<<",x="<<x<<",y="<<y<<",xprim="<<xprime[counter]<<endl;

        counter++;
    }

    double minval=1000;
    int minIndex=-1;

    //find the leftmost point
    for (int i=1; i<counter; i++){
        if (xprime[i]<minval)
        {
            minval=xprime[i];
            minIndex=i;
        }
    }

    cout<<"minval= "<<minval<<endl;

    double xm_new[2];

    if (minIndex==-1){
        xm_new[0]=xc2;
        xm_new[1]=yc2;
    }else {
        xm_new[0] = xc2 + 0.8 * (a + minval) * cos(theta2);
        xm_new[1] = yc2 + 0.8 * (a + minval) * sin(theta2);
    }

    for (int i=1; i<counter; i++){
        tempAd[i][0]=tempAd[i][0]+xc2-xm_new[0];
        tempAd[i][1]=tempAd[i][1]+yc2-xm_new[1];
    }

    for (int i=1; i<=totalnum2; i++){

        if (i<counter){
            Adhesion2[i][0]=tempAd[i][0];
            Adhesion2[i][1]=tempAd[i][1];
            AdhesionIndex2[i][0]=tempAdIndex[i][0];
            AdhesionIndex2[i][1]=tempAdIndex[i][1];
            Adhesion_shift2[i][0]=tempAdshift[i][0];
            Adhesion_shift2[i][1]=tempAdshift[i][1];
        }else{
            Adhesion2[i][0]=0.0;
            Adhesion2[i][1]=0.0;
            AdhesionIndex2[i][0]=0;
            AdhesionIndex2[i][1]=0;
            Adhesion_shift2[i][0]=0;
            Adhesion_shift2[i][1]=0;
        }

        ifremove2[i]=0;
        Adhesion02[i][0]=0;
        Adhesion02[i][1]=0;
    }

    xm2[0]=xm_new[0];
    xm2[1]=xm_new[1];

    //reset the totalnum
    totalnum2=counter-1;

    //free the tempertory pointer
    delete[] tempAd;

    delete[] tempAdIndex;

    delete[] tempAdshift;
    //--------------------------
}
*/
