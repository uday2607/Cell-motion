//
// Created by jingchen Feng on 4/11/17.
//

bool ifcell(double x, double y, double xm, double ym, double a, double b, double theta) {
    double val;
    bool flag;

    double xcan, ycan;

    xcan = (x - xm) * cos(theta) + (y - ym) * sin(theta);
    ycan = -(x - xm) * sin(theta) + (y - ym) * cos(theta);

    val = xcan * xcan / (a * a) + ycan * ycan / (b * b);

    flag = (val <= 1);

    return flag;
}

bool
ifcell_periodic(double x, double y, double xm, double ym, double a, double b, double theta, double shift[], int direc) {
    //judge if a node (rx, ry) is inside the cell, given periodic boundary condition

    double x0 = x;
    double y0 = y;

    bool isin = false;

    if (ifcell(x, y, xm, ym, a, b, theta)) {
        isin = true;
        shift[0] = 0.0;
        shift[1] = 0.0;
        return isin;
    } else {
        //depending on which boundary cell is moving out, we determine the direction of shift
        if (direc == 0) {
            isin = false;
            shift[0] = 0.0;
            shift[1] = 0.0;   //this node is not inside the cell
            return isin;
        }

        if (direc == 1) {

            //case 1, x=x0+Lx
            if (ifcell(x0 + Lx, y, xm, ym, a, b, theta)) {
                isin = true;
                shift[0] = Lx;
                shift[1] = 0.0;
                return isin;
            }

            //case 2, x=x0+Lx, y=y0+Ly*sqrt(3)/2  x=x0+1/2*Ly
            if (ifcell(x0 + Lx + Ly / 2, y0 + h, xm, ym, a, b, theta)) {
                isin = true;
                shift[0] = Lx + Ly / 2;
                shift[1] = h;  //h=Ly*0.86602540378443864676*(1+eps);
                return isin;
            }

            //case 3, x=x0+Lx, y=y0-Ly*sqrt(3)/2 x=x0-1/2*Ly
            if (ifcell(x0 + Lx - Ly / 2, y0 - h, xm, ym, a, b, theta)) {
                isin = true;
                shift[0] = Lx - Ly / 2;
                shift[1] = -h;
                return isin;
            }

            //case 4, x=x0-Lx
            if (ifcell(x0 - Lx, y, xm, ym, a, b, theta)) {
                isin = true;
                shift[0] = -Lx;
                shift[1] = 0.0;
                return isin;
            }

            //case 5, x=x0-Lx, y=y0+Ly*sqrt(3)/2 x=x0+1/2*Ly
            if (ifcell(x0 - Lx + Ly / 2, y0 + h, xm, ym, a, b, theta)) {
                isin = true;
                shift[0] = -Lx + Ly / 2;
                shift[1] = h;
                return isin;
            }

            //case 6, x=x0-Lx, y=y0-Ly*sqrt(3)/2 x=x0-1/2*Ly
            if (ifcell(x0 - Lx - Ly / 2, y0 - h, xm, ym, a, b, theta)) {
                isin = true;
                shift[0] = -Lx - Ly / 2;
                shift[1] = -h;
                return isin;
            }

            //case 7, y=y0+Ly*sqrt(3)/2 x=x0+1/2*Ly
            if (ifcell(x0 + Ly / 2, y0 + h, xm, ym, a, b, theta)) {
                isin = true;
                shift[0] = Ly / 2;
                shift[1] = h;
                return isin;
            }

            //case 8, y=y0-Ly*sqrt(3)/2 x=x0-1/2*Ly
            if (ifcell(x0 - Ly / 2, y0 - h, xm, ym, a, b, theta)) {
                isin = true;
                shift[0] = -Ly / 2;
                shift[1] = -h;
                return isin;
            }


        }

    }

    shift[0] = 0.0;
    shift[1] = 0.0;
    return isin;
}