function retval = wertex( si, sj, mindist )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculates pseudo vertex from one contact and two Pmax spheres
%
% si: Pmax sphere cell i
% sj: Pmax sphere cell j
%
% each sphere organized as
% s(1) = R
% s(2) = Mx
% s(3) = My
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find intersection points
Ri = si(1);
Rj = sj(1);
Mi = si(2:3);
Mj = sj(2:3);
pts = intsph( Mi,Ri, Mj,Rj, mindist );

sizpts = size(pts);
switch( sizpts(1) )
  case 3
    error('stop wertex.m: spheres identical');
  case 2
    retval = pts;
  case 1
    retval = pts;
  case 0
    warning('wertex.m: spheres not intersecting');
end

