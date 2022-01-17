# Source: https://link.springer.com/article/10.1007/s00791-013-0214-3
# From: https://github.com/chraibi/EEOver/blob/master/solvers.c
include("Poly_solvers.jl")

function ell2tr(x, y, aa, bb, cc,
        dd, ee, ff)

    return aa*x*x + bb*x*y + cc*y*y + dd*x + ee*y + ff    
end

function nointpts(a1, b1, a2, b2, h1,
    k1, h2, k2, t1, t2, h2_tr,
    k2_tr, aa, bb, cc, dd, ee,
    ff)

    # some tmp-variables to avoid doing things several times    
    a1b1 = a1*b1 
    a2b2 = a2*b2 
    area_1 = pi*a1b1
    area_2 = pi*a2b2 
    EPS = 10.0^-5

    # The relative size of the two ellipses can be found from the axis
    # lengths 
    relsize = a1b1 - a2b2 
    if (relsize > 0.0)
        # First Ellipse is larger than second ellipse.
        # If second ellipse center (H2_TR, K2_TR) is inside
        # first ellipse, then ellipse 2 is completely inside 
        # ellipse 1. Otherwise, the ellipses are disjoint.

        if ((h2_tr*h2_tr)/(a1*a1) + 
            (k2_tr*k2_tr)/(b1*b1) < 1.0)
            return area_2
        else
            return 0.0    
        end    

    elseif (relsize < 0.0)    
        # Second Ellipse is larger than first ellipse
        # If first ellipse center (0, 0) is inside the
        # second ellipse, then ellipse 1 is completely inside
        # ellipse 2. Otherwise, the ellipses are disjoint
        # AA*x^2 + BB*x*y + CC*y^2 + DD*x + EE*y + FF = 0

        if (ff < 0.0)
            return area_1 
        else 
            return 0.0
        end
    
    else 
        # If execution arrives here, the relative sizes are identical.
        # Are the ellipses the same?  Check the parameters to see.
        # MC. Ellipses are the same if: H1=H2 And K1==K2 And Area_1 == Area_2
        if (abs(h1 - h2) < EPS && abs(k1 - k2) < EPS && abs(area_1 - area_2) < EPS)
            return area_1 
        else     
            return 0.0
        end     
    end    
end    

# two distinct intersection points (x1, y1) and (x2, y2) find overlap area
function twointpts(x, y, a1, b1, t1,
    a2, b2, h2_tr, k2_tr, t2, aa,
    bb, cc, dd, ee, ff)

    #temporary constant 
    EPS = 10.0^-5
    
    # if execution arrives here, the intersection points are not
    # tangents.
    # determine which direction to integrate in the ellipse_segment
    # routine for each ellipse.
    # find the parametric angles for each point on ellipse 1
    if (abs(x[1]) > a1)
        if x[1] < 0
           x[1] = -a1
        else
           x[1] = a1
        end        
    end
    
    if (y[1] < 0.0) 
        # Quad III or IV
        theta1 = 2.0*pi - acos(x[1]/a1)
    else
        # Quad I or II
        theta1 = acos(x[1]/a1)
    end

    if (abs(x[2]) > a1)
        if x[2] < 0
           x[2] = -a1
        else
           x[2] = a1
        end        
    end
    
    if (y[2] < 0.0) 
        # Quad III or IV
        theta2 = 2.0*pi - acos(x[2]/a1)
    else
        # Quad I or II
        theta2 = acos(x[2]/a1)
    end

    # logic is for proceeding counterclockwise from theta1 to theta2
    if theta1 > theta2 
        tmp = theta1 
        theta1 = theta2 
        theta2 = tmp
    end 
    
    # find a point on the first ellipse that is different than the two
    # intersection points.
    xmid = a1*cos((theta1 + theta2)/2.0)
    ymid = b1*sin((theta1 + theta2)/2.0)

    # the point (xmid, ymid) is on the first ellipse 'between' the two
    # intersection points (x[1], y[1]) and (x[2], y[2]) when travelling 
    # counter- clockwise from (x[1], y[1]) to (x[2], y[2]).  If the point
    # (xmid, ymid) is inside the second ellipse, then the desired segment
    # of ellipse 1 contains the point (xmid, ymid), so integrate 
    # counterclockwise from (x[1], y[1]) to (x[2], y[2]).  Otherwise, 
    # integrate counterclockwise from (x[2], y[2]) to (x[1], y[1])
    if (ell2tr(xmid, ymid, aa, bb, cc, dd, ee, ff) > 0.0)
        tmp = theta1 
        theta1 = theta2 
        theta2 = tmp
    end    

    # here is the ellipse segment routine for the first ellipse
    if (theta1 > theta2)
        theta1 -= 2.0*pi 
    end 
    if ((theta2 - theta1) > pi)
        trsign = 1.0
    else 
        trsign = -1.0
    end 
    
    area1 = 0.5*(a1*b1*(theta2-theta1) + trsign*abs(x[1]y[2] - x[2]y[1]))

    if (area1 < 0)
        area1 += a1*b1
    end
    
    # find ellipse 2 segment area.  The ellipse segment routine
    # needs an ellipse that is centered at the origin and oriented
    # with the coordinate axes.  The intersection points (x[1], y[1]) and
    # (x[2], y[2]) are found with both ellipses translated and rotated by
    # (-H1, -K1) and -PHI_1.  Further translate and rotate the points
    # to put the second ellipse at the origin and oriented with the
    # coordinate axes.  The translation is (-H2_TR, -K2_TR), and the
    # rotation is -(PHI_2 - PHI_1) = PHI_1 - PHI_2
    cosphi = cos(t1 - t2)
    sinphi = sin(t1 - t2)
    x1_tr = (x[1] - h2_tr)*cosphi + (y[1] - k2_tr)*(-1.0*sinphi)
    y1_tr = (x[1] - h2_tr)*sinphi + (y[1] - k2_tr)*cosphi 
    x2_tr = (x[2] - h2_tr)*cosphi + (y[2] - k2_tr)*(-1.0*sinphi)
    y2_tr = (x[2] - h2_tr)*sinphi + (y[2] - k2_tr)*cosphi

    # determine which branch of the ellipse to integrate by finding a
    # point on the second ellipse, and asking whether it is inside the
    # first ellipse (in their once-translated+rotated positions)
    # find the parametric angles for each point on ellipse 1
    if (abs(x1_tr) > a2)
        if x1_tr < 0
           x1_tr = -a2
        else
           x1_tr = a2
        end        
    end
    
    if (y1_tr < 0.0) 
        # Quad III or IV
        theta1 = 2.0*pi - acos(x1_tr/a2)
    else
        # Quad I or II
        theta1 = acos(x1_tr/a2)
    end

    if (abs(x2_tr) > a2)
        if x2_tr < 0
           x2_tr = -a2
        else
           x2_tr = a2
        end        
    end
    
    if (y2_tr < 0.0) 
        # Quad III or IV
        theta2 = 2.0*pi - acos(x2_tr/a2)
    else
        # Quad I or II
        theta2 = acos(x2_tr/a2)
    end  

    # logic is for proceeding counterclockwise from theta1 to theta2
    if theta1 > theta2 
        tmp = theta1 
        theta1 = theta2 
        theta2 = tmp
    end 

    # find a point on the second ellipse that is different than the two
    # intersection points.
    xmid = a2*cos((theta1+theta2)/2.0)
    ymid = b2*sin((theta1+theta2)/2.0)

    # translate the point back to the second ellipse in its once-
    # translated+rotated position
    cosphi = cos(t2 - t1)
    sinphi = sin(t2 - t1)
    xmid_rt = xmid*cosphi + ymid*(-1.0*sinphi) + h2_tr
    ymid_rt = xmid*sinphi + ymid*cosphi + k2_tr 

    # the point (xmid_rt, ymid_rt) is on the second ellipse 'between' the
    # intersection points (x[1], y[1]) and (x[2], y[2]) when travelling
    # counterclockwise from (x[1], y[1]) to (x[2], y[2]).  If the point
    # (xmid_rt, ymid_rt) is inside the first ellipse, then the desired 
    # segment of ellipse 2 contains the point (xmid_rt, ymid_rt), so 
    # integrate counterclockwise from (x[1], y[1]) to (x[2], y[2]).  
    # Otherwise, integrate counterclockwise from (x[2], y[2]) to 
    # (x[1], y[1])
    if (((xmid_rt*xmid_rt)/(a1*a1) + (ymid_rt*ymid_rt)/(b1*b1)) > 1.0)
        tmp = theta1 
        theta1 = theta2 
        theta2 = tmp 
    end

    # here is the ellipse segment routine for the second ellipse
    if (theta1 > theta2)
        theta1 -= 2.0*pi 
    end 
    if ((theta2 - theta1) > pi)
        trsign = 1.0
    else 
        trsign = -1.0
    end 

    area2 = 0.5*(a2*b2*(theta2 - theta1) +  trsign*abs(x1_tr*y2_tr - x2_tr*y1_tr))

    if (area2 < 0)
        area2 += a2*b2 
    end
    
    return area1 +  area2

end    

# check whether an intersection point is a tangent or a cross-point
function istanpt(x, y, a1, b1, aa, 
        bb, cc, dd, ee, ff)

    # Avoid inverse trig calculation errors: there could be an error 
    # if |x1/A| > 1.0 when calling acos().  If execution arrives here, 
    # then the point is on the ellipse within EPS.   
    if (abs(x) > a1)
        if x < 0
            x = -1.0*a1 
        else
            x = a1 
        end        
    end
    
    # Calculate the parametric angle on the ellipse for (x, y)
    # The parametric angles depend on the quadrant where each point
    # is located.
    if (y < 0.0)
        theta = 2.0*pi - acos(x/a1)
    else 
        theta = acos(x/a1)
    end 
    
    eps_radian = 0.1 #arbitrary value 

    # determine two points that are on each side of (x, y) and lie on
    # the first ellipse
    x1 = a1*cos(theta + eps_radian)
    y1 = b1*sin(theta + eps_radian)
    x2 = a1*cos(theta - eps_radian)
    y2 = b1*sin(theta - eps_radian)

    # evaluate the two adjacent points in the second ellipse equation
    test1 = ell2tr(x1, y1, aa, bb, cc, dd, ee, ff)
    test2 = ell2tr(x2, y2, aa, bb, cc, dd, ee, ff)

    # if the ellipses are tangent at the intersection point, then
    # points on both sides will either both be inside ellipse 1, or
    # they will both be outside ellipse 1

    if ((test1*test2) > 0.0)
        return 1
    else 
        return 0
    end         
end         

function threeintpts(xint, yint, a1, b1, t1,
    a2, b2, h2_tr, k2_tr, t2, aa,
    bb, cc, dd, ee, ff)

    # need to determine which point is a tangent, and which two points
    # are intersections
    tanpts = 0
    tanindex = 0
    for i in 1:3
        fnRtn = istanpt(xint[i], yint[i], a1, b1, aa, bb, cc, dd, ee, ff)

        if (fnRtn == 1)
            tanpts += 1
            tanindex = i
        end
    end    

    # there MUST be 2 intersection points and only one tangent
    if (tanpts != 1)
        # some weird issue. should never even happen
        return -1.0
    end    

    # store the two interesection points into (x[1], y[1]) and 
    # (x[2], y[2])
    if (tanindex == 0)
        xint[1] = xint[3]
        yint[1] = yint[3]
    elseif tanindex == 1
        xint[2] = xint[3]
        yint[2] = yint[3]
    end    

    area = twointpts(xint, yint, a1, b1, t1, a2, b2, h2_tr, k2_tr,
                t2, aa, bb, cc, dd, ee, ff)
    
    return area            

end    

# four intersection points
function fourintpts(xint, yint, a1, b1, t1,
    a2, b2, h2_tr, k2_tr, t2, aa,
    bb, cc, dd, ee, ff)

    # some tmp-variables to avoid calculating the same thing several times
    a1b1 = a1*b1 
    a2b2 = a2*b2
    area_1 = pi*a1b1
    area_2 = pi*a2b2 

    # only one case, which involves two segments from each ellipse, plus
    # two triangles.
    # get the parametric angles along the first ellipse for each of the
    # intersection points

    theta = []
    for i in 1:4
        if abs(xint[i]) > a1 
            if xint[i] < 0
                xint[i] = -1.0*A1 
            else 
                xint[i] = a1 
            end        
        end
        
        if (yint[i] < 0.0)
            push!(theta, 2.0*pi - acos(xint[i]/a1))
        else
            push!(theta, acos(xint[i]/a1))  
        end       
    end
    
    # sort the angles by straight insertion, and put the points in 
    # counter-clockwise order
    for j in 1:3
        tmp0 = theta[j]
        tmp1 = xint[j]
        tmp2 = yint[j]

        for k in reverse(1:j)
            
            if (theta[k] <= tmp0)
                break
            end
            
            theta[k+1] = theta[k]
            xint[k+1] = xint[k]
            yint[k+1] = yint[k]
        end
        
        theta[1] = tmp0
        xint[1] = tmp1
        yint[1] = tmp2 
    end

    # find the area of the interior quadrilateral
    area1 = 0.5*abs((xint[3] - xint[1])*(yint[4] - yint[2]) - 
                    (xint[4] - xint[2])*(yint[3] - yint[1]))

    # the intersection points lie on the second ellipse in its once
    # translated+rotated position.  The segment algorithm is implemented
    # for an ellipse that is centered at the origin, and oriented with
    # the coordinate axes; so, in order to use the segment algorithm
    # with the second ellipse, the intersection points must be further
    # translated+rotated by amounts that put the second ellipse centered
    # at the origin and oriented with the coordinate axes.
    cosphi = cos(t1 - t2)
    sinphi = sin(t1 - t2)

    xint_tr = Real[0.0, 0.0, 0.0, 0.0]
    yint_tr = Real[0.0, 0.0, 0.0, 0.0]
    theta_tr = Real[0.0, 0.0, 0.0, 0.0]
    for i in 1:4

        xint_tr[i] = (xint[i] - h2_tr)*cosphi + (yint[i] - k2_tr)*(-1.0*sinphi)
        yint_tr[i] = (xint[i] - h2_tr)*sinphi + (yint[i] - k2_tr)*cosphi 

        if (abs(xint_tr[i]) > a2)
            if xint_tr[i] < 0
                xint_tr[i] = -1.0*a2 
            else
                xint_tr = a2    
            end    
        end
        
        if (yint_tr[i] < 0.0)
            theta_tr[i] = 2.0*pi - acos(xint_tr[i]/a2)
        else 
            theta_tr[i] = acos(xint_tr[i]/a2)
        end 
    end    
    
    # get the area of the two segments on ellipse 1
    xmid = a1*cos((theta[1] + theta[2])/2.0)
    ymid = b1*sin((theta[1] + theta[2])/2.0)  

    # the point (xmid, ymid) is on the first ellipse 'between' the two
    # sorted intersection points (xint[1], yint[1]) and (xint[2], yint[2])
    # when travelling counter- clockwise from (xint[1], yint[1]) to 
    # (xint[2], yint[2]).  If the point (xmid, ymid) is inside the second 
    # ellipse, then one desired segment of ellipse 1 contains the point 
    # (xmid, ymid), so integrate counterclockwise from (xint[1], yint[1])
    # to (xint[2], yint[2]) for the first segment, and from 
    # (xint[3], yint[3] to (xint[4], yint[4]) for the second segment.   
    if (ell2tr(xmid, ymid, aa, bb, cc, dd, ee, ff) < 0.0)
        area2 = 0.5*(a1b1*(theta[2]-theta[1]) - abs(xint[1]*yint[2] - xint[2]*yint[1]))
        area3 = 0.5*(a1b1*(theta[4]-theta[3]) - abs(xint[3]*yint[4] - xint[4]*yint[3]))
        area4 = 0.5*(a2b2*(theta_tr[3] - theta_tr[2]) - abs(xint_tr[2]*yint_tr[3] - xint_tr[3]*yint_tr[2]))

        if (theta_tr[4] > theta_tr[1])
            area5 = 0.5*(a2b2*(theta_tr[1] - (theta_tr[4] - 2.0*pi))
                    - abs(xint_tr[4]*yint_tr[1] - xint_tr[1]*yint_tr[4]))
        else 
            area5 = 0.5*(a2b2*(theta_tr[1] - theta_tr[4])
                    - abs(xint_tr[4]*yint_tr[1] - xint_tr[1]*yint_tr[4]))
        end         
    
    else
        area2 = 0.5*(a1b1*(theta[3]-theta[2]) - abs(xint[2]*yint[3] - xint[3]*yint[2]))
        area3 = 0.5*(a1b1*(theta[1]-(theta[4]-2.0*pi)) - abs(xint[4]*yint[1] - xint[1]*yint[4]))
        area4 = 0.5*(a2b2*(theta_tr[2] - theta_tr[1]) - abs(xint_tr[1]*yint_tr[2] - xint_tr[2]*yint_tr[1]))
        area5 = 0.5*(a2b2*(theta_tr[4] - theta_tr[3])
                    - abs(xint_tr[3]*yint_tr[4] - xint_tr[4]*yint_tr[3]))
    end

    if (area5 < 0.0)
        area5 += area_2 
    end 
    
    if (area4 < 0.0)
        area4 += area_2 
    end 

    if (area3 < 0.0)
        area3 += area_1 
    end 

    if (area2 < 0.0)
        area2 += area_1 
    end 

    # Overlap area 
    Area = area1 + area2 + area3 + area4 + area5
    return Area

end    

function ellipse_ellipse_overlap(a1, b1, h1,
                    k1, t1, a2, b2,
                    h2, k2, t2)
    
    # Useful Constants
    twopi = 2.0*pi
    EPS = 10.0^-5

    # Angles should be in the range (-pi, pi)
    if abs(t1) > pi
        t1 = mod(t1, pi)
    end 
    if abs(t2) > pi 
        t2 = mod(t2, pi)
    end 
    
    # ==================================================================
    # = DETERMINE THE TWO ELLIPSE EQUATIONS FROM INPUT PARAMETERS =
    # ===================================================================
    # Finding the points of intersection between two general ellipses
    # requires solving a quartic equation.  Before attempting to solve the
    # quartic, several quick tests can be used to eliminate some cases
    # where the ellipses do not intersect.  Optionally, can whittle away
    # at the problem, by addressing the easiest cases first.
    # Working with the translated+rotated ellipses simplifies the
    # calculations.  The ellipses are translated then rotated so that the
    # first ellipse is centered at the origin and oriented with the 
    # coordinate axes.  Then, the first ellipse will have the implicit 
    # (polynomial) form of
    #   x^2/A1^2 + y+2/B1^2 = 1

    # For the second ellipse, the center is first translated by the amount
    # required to put the first ellipse at the origin, e.g., by (-H1, -K1)  
    # Then, the center of the second ellipse is rotated by the amount
    # required to orient the first ellipse with the coordinate axes, e.g.,
    # through the angle -PHI_1.
    # The translated and rotated center point coordinates for the second
    # ellipse are found with the rotation matrix, derivations are 
    # described in the reference.
    cosphi = cos(t1)
    sinphi = sin(t1)
    h2_tr = (h2 - h1)*cosphi + (k2-k1)*sinphi 
    k2_tr = (h1 - h2)*sinphi + (k2-k1)*cosphi
    phi = t2 - t1 
    if abs(phi) > twopi
        phi = mod(phi, twopi)
    end

    # Calculate implicit (Polynomial) coefficients for the second ellipse
    # in its translated-by (-H1, -H2) and rotated-by -PHI_1 postion
    #       AA*x^2 + BB*x*y + CC*y^2 + DD*x + EE*y + FF = 0
    # Formulas derived in the reference
    # To speed things up, store multiply-used expressions first
    cosphi = cos(t2)
    cosphi2 = cosphi*cosphi 
    sinphi = sin(t2)
    sinphi2 = sinphi*sinphi
    cosphisinphi = 2.0*cosphi*sinphi 
    aa2 = a2*a2
    bb2 = b2*b2 
    tmp0 = (cosphi*h2_tr + sinphi*k2_tr)/aa2
    tmp1 = (sinphi*h2_tr - cosphi*k2_tr)/bb2
    tmp2 = cosphi*h2_tr + sinphi*k2_tr
    tmp3 = sinphi*h2_tr - cosphi*k2_tr

    # implicit polynomial coefficients for the second ellipse
    aa = cosphi2/aa2 + sinphi2/bb2 
    bb = cosphisinphi/aa2 - cosphisinphi/bb2 
    cc = sinphi2/aa2 + cosphi2/bb2
    dd = -2.0*cosphi*tmp0 - 2.0*sinphi*tmp1
    ee = -2.0*sinphi*tmp0 + 2.0*cosphi*tmp1
    ff = tmp2*tmp2/aa2 + tmp3*tmp3/bb2 - 1.0

    # =======================================================================
    # = CREATE AND SOLVE THE QUARTIC EQUATION TO FIND INTERSECTION POINTS ==
    # =======================================================================
    # If execution arrives here, the ellipses are at least 'close' to
    # intersecting.
    # Coefficients for the Quartic Polynomial in y are calculated from
    # the two implicit equations.
    # Formulas for these coefficients are derived in the reference.
    c_5 = a1*a1*a1*a1*aa*aa + b1*b1*(a1*a1*(bb*bb - 2.0*aa*cc) + b1*b1*cc*cc)
    c_4 = 2.0*b1*(b1*b1*cc*ee + a1*a1*(bb*dd - aa*ee))
    c_3 = (a1*a1*((b1*b1*(2.0*aa*cc - bb*bb) + dd*dd - 2.0*aa*ff) - 2.0*a1*a1*aa*aa) + 
          b1*b1*(2.0*cc*ff + ee*ee))
    c_2 = 2.0*b1*(a1*a1*(aa*ee - bb*dd) + ee*ff)
    c_1 = (a1*(a1*aa - dd) + ff)*(a1*(a1*aa + dd) + ff)
    
    # Once the coefficients for the Quartic Equation in y are known, the
    # roots of the quartic polynomial will represent y-values of the 
    # intersection points of the two ellipse curves.
    # The quartic sometimes degenerates into a polynomial of lesser 
    # degree, so handle all possible cases.

    # Depending on the coeffs, find the roots (only real)
    if c_5 != 0.0
        poly_roots = ferrari(c_5, c_4, c_3, c_2, c_1)
    elseif c_4 != 0.0
        poly_roots = cardan(c_4, c_3, c_2, c_2)
    elseif c_3 != 0.0
        poly_roots = roots2(c_3, c_2, c_1)
    elseif c_2 != 0.0
        poly_roots = ComplexF64[-1.0*c_1/c_2]  
    else
        poly_roots = ComplexF64[]              
    end    

    # Get the real roots 
    ychk = []
    for rt in poly_roots
        if isreal(rt)
            push!(ychk, Real(rt*b1))
        end
    end
    
    # reverse sort the roots 
    ychk = sort(ychk, rev=true)
    nychk = size(ychk)[1]

    # determine whether the polynomial roots are points of 
    # intersection for the two ellipses 
    nintpts = 0
    xint = []
    yint = []
    for i in 1:nychk

        # check for multiple roots 
        if (i < nychk && abs(ychk[i] - ychk[i+1]) < EPS)
            continue 
        end 
        
        # check intersection points for ychk[i]
        if (abs(ychk[i]) > b1)
            x1 = 0.0
        else 
            x1 = a1*sqrt(1.0 - (ychk[i]*ychk[i])/(b1*b1)) 
        end 
        x2 = -1.0*x1
        
        if (abs(ell2tr(x1, ychk[i], aa, bb, cc, dd, ee, ff)) < EPS)
            nintpts += 1

            if (nintpts > 4)
                return -1.0      
            end
            push!(xint, x1) 
            push!(yint, ychk[i])
        end            

        if (abs(ell2tr(x2, ychk[i], aa, bb, cc, dd, ee, ff)) < EPS
            && abs(x2 - x1) > EPS)
            nintpts += 1

            if (nintpts > 4)
                return -1.0      
            end
            push!(xint, x2) 
            push!(yint, ychk[i])
        end
    end

    # ====================================================================
    # = HANDLE ALL CASES FOR THE NUMBER OF INTERSCTION POINTS ============
    # ========================================xs==========================
    if nintpts == 1 || nintpts == 0
        Area = nointpts(a1, b1, a2, b2, h1, k1, h2, k2, t1, t2, 
                h2_tr, k2_tr, aa, bb, cc, dd, ee, ff)       
        
        return Area
    elseif nintpts == 2    

        # when there are two intersection points, it is possible for
        # them to both be tangents, in which case one of the ellipses
        # is fully contained within the other.  Check the points for
        # tangents; if one of the points is a tangent, then the other
        # must be as well, otherwise there would be more than 2 
        # intersection points.     
        fnRtnCode = istanpt(xint[1], yint[1], a1, b1, aa, bb, cc, dd, ee, ff)
        
        if (fnRtnCode == 1)
            Area = nointpts(a1, b1, a2, b2, h1, k1, h2, k2, t1, t2, 
            h2_tr, k2_tr, aa, bb, cc, dd, ee, ff)
        else
            Area = twointpts(xint, yint, a1, b1, t1, a2, b2, h2_tr, 
                    k2_tr, t2, aa, bb, cc, dd, ee, ff)    
        end
        
        return Area 

    elseif nintpts == 3
        # when there are three intersection points, one and only one
        # of the points must be a tangent point.
        Area = threeintpts(xint, yint, a1, b1, t1, a2, b2, h2_tr,
                k2_tr, t2, aa, bb, cc, dd, ee, ff)
        return Area 
    
    elseif nintpts == 4
        # four intersections points has only one case.
        Area = fourintpts(xint, yint, a1, b1, t1, a2, b2, h2_tr,
                k2_tr, t2, aa, bb, cc, dd, ee, ff)    
        return Area        
    else
        println(nintpts)
        return -1.0    
    end
end