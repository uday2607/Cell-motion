function retval = arc( R, M, angran, numpts, minang )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% arc.m
%%%%%%%%
%%%%%%%% returns a matrix of points on circular arc
%%%%%%%%
%%%%%%%% R                             circle radius
%%%%%%%% M = [ Mx, My ]                circle center
%%%%%%%% angran = [ thmin, thmax ]     range of angles (within -pi,pi)
%%%%%%%% numpts                        number of points
%%%%%%%% 
%%%%%%%% Care is taken that the correct angle range is returned.
%%%%%%%% In particular, once thmin > 0, and thmax < 0, 2*pi is 
%%%%%%%% added to thmax in order to obtain the appropriate arc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thmin = angran(1);
if angran(1)>0 && angran(2)<0
  thmax = angran(2) + 2*pi;
else
  thmax = angran(2);
end

if thmax-thmin < minang
  warning('arc.m: very short arc, approximating by single point');
  theta = thmin;
else
  theta = thmin:(thmax-thmin)/(numpts-1):thmax;
end

x = M(1) + R*cos(theta);
y = M(2) + R*sin(theta);

retval = [ x; y ];

