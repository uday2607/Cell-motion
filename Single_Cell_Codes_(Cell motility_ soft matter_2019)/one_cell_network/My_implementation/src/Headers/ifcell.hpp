#include <cmath>

bool if_node_in_cell(double x, double y, double xm, double ym,
                     double a, double b, double theta) {
    // check if a node is in the cell

    double val, xcan, ycan;
    bool flag;

    // Cell coordinate from the center of the cell
    xcan = (x - xm)*cos(theta) + (y - ym)*sin(theta);
    ycan = -(x - xm)*sin(theta) + (y - ym)*cos(theta);

    // Ellipse equation
    val = xcan * xcan / (a * a) + ycan * ycan / (b * b);

    // Check if the node is inside the cell
    flag = (val <= 1);

    return flag;
}

bool if_cell_is_periodic(double x, double y, double xm,
  double ym, double a, double b, double theta, double shift[], int direc) {
    // check if a node (rx, ry) is inside the cell,
    // given periodic boundary condition

    double x0 = x;
    double y0 = y;

    // flag to check
    bool isin = false;

    // If the node is in the cell, don't shift
    if (if_node_in_cell(x, y, xm, ym, a, b, theta)) {

        isin = true;
        shift[0] = 0.0;
        shift[1] = 0.0;

        return isin;

    } else {
        // depending on which boundary cell is moving out,
        // we determine the direction of shift
        if (direc == 0) {

            isin = false;
            shift[0] = 0.0;
            shift[1] = 0.0;   //this node is not inside the cell

            return isin;

        }

        if (direc == 1) {

            //case 1, x=x0+Lx
            if (if_node_in_cell(x0 + Lx, y, xm, ym, a, b, theta)) {

                isin = true;
                shift[0] = Lx;
                shift[1] = 0.0;

                return isin;

            }

            //case 2, x=x0+Lx, y=y0+Ly*sqrt(3)/2  x=x0+1/2*Ly
            if (if_node_in_cell(x0 + Lx + Ly / 2, y0 + h, xm,
                                ym, a, b, theta)) {

                isin = true;
                shift[0] = Lx + Ly / 2;
                shift[1] = h;  //h=Ly*0.86602540378443864676*(1+eps);

                return isin;
            }

            //case 3, x=x0+Lx, y=y0-Ly*sqrt(3)/2 x=x0-1/2*Ly
            if (if_node_in_cell(x0 + Lx - Ly / 2, y0 - h, xm, ym, a, b, theta)) {

                isin = true;
                shift[0] = Lx - Ly / 2;
                shift[1] = -h;

                return isin;
            }

            //case 4, x=x0-Lx
            if (if_node_in_cell(x0 - Lx, y, xm, ym, a, b, theta)) {

                isin = true;
                shift[0] = -Lx;
                shift[1] = 0.0;

                return isin;
            }

            //case 5, x=x0-Lx, y=y0+Ly*sqrt(3)/2 x=x0+1/2*Ly
            if (if_node_in_cell(x0 - Lx + Ly / 2, y0 + h, xm, ym, a, b, theta)) {

                isin = true;
                shift[0] = -Lx + Ly / 2;
                shift[1] = h;

                return isin;
            }

            //case 6, x=x0-Lx, y=y0-Ly*sqrt(3)/2 x=x0-1/2*Ly
            if (if_node_in_cell(x0 - Lx - Ly / 2, y0 - h, xm, ym, a, b, theta)) {

                isin = true;
                shift[0] = -Lx - Ly / 2;
                shift[1] = -h;

                return isin;
            }

            //case 7, y=y0+Ly*sqrt(3)/2 x=x0+1/2*Ly
            if (if_node_in_cell(x0 + Ly / 2, y0 + h, xm, ym, a, b, theta)) {

                isin = true;
                shift[0] = Ly / 2;
                shift[1] = h;

                return isin;
            }

            //case 8, y=y0-Ly*sqrt(3)/2 x=x0-1/2*Ly
            if (if_node_in_cell(x0 - Ly / 2, y0 - h, xm, ym, a, b, theta)) {

                isin = true;
                shift[0] = -Ly / 2;
                shift[1] = -h;

                return isin;
            }


        }

    }
    // If none of the cases are true...
    shift[0] = 0.0;
    shift[1] = 0.0;
    return isin;
}
