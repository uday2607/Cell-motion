function retval = intsph( C1,R1, C2,R2, mindist )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% intsph.m
%%%%%%%%
%%%%%%%% intersect two spheres 2d
%%%%%%%% C1 center sphere 1 (row vector [x,y]), R1 radius sphere 1
%%%%%%%% C2 center sphere 2 (row vector [x,y]), R2 radius sphere 2
%%%%%%%% mindist: minimal distance (when to consider two points equal)
%%%%%%%%
%%%%%%%% Remark:
%%%%%%%% spheres intersecting iff
%%%%%%%% |R1-R2| <= |C1-C2| <= R1+R2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

deltR = R1 - R2;
summR = R1 + R2;
deltC = C1 - C2;
dis2C = deltC*deltC';
distC = sqrt( dis2C );

% spheres not intersecting
if abs(deltR) > distC || summR < distC
  retval = [];
  return;
end

% spheres coinciding
if abs(deltR) < mindist && distC < mindist
  retval = ones(3,2);    % three pts represent infinitely many
  return;
end

% direct local unit vector towards bigger sphere
if R1 > R2
  B = R1;
  A = R2;
  bc = C1;
  ac = C2;
else
  B = R2;
  A = R1;
  bc = C2;
  ac = C1;
end
deltC =   bc - ac;
orig  = ( ac + bc )/2;

% spheres barely in contact
if abs( distC - summR ) < mindist || abs( distC - abs(deltR) ) < mindist
  y = 0;
  x = distC/2 - B;
else
% spheres normal
  xx = ( A^2 - B^2 )/( 2*distC );
  pb = distC/2 - xx;
  yy = sqrt( B^2 - pb^2 );
  ym = -yy;
  yp =  yy;
  x = [ xx; xx ];
  y = [ ym; yp ];
end

Dphi = atan2( deltC(2), deltC(1) );
xy = ShiftRot2d( x,y, orig(1),orig(2), Dphi );

retval = xy;

