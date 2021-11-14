function retval = vertex( s1, s2, s3, c1, c2, c3, mindist )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculates vertices from three contact spheres
% - usually one vertex point is returned.
% - in special colinear configurations it might be that both vertex
%   candidates are equally far away from the cells.
%   then both are returned.
% - when three cells colinear another possibility is, that there is
%   no vertex at all. then an empty vector is returned.
%
%
% s1: sphere1
% s2: sphere2
% s3: sphere3
% c1: cell1
% c2: cell2
% c3: cell3
%
% each sphere organized as
% s(1) = Rij
% s(2) = Mijx
% s(3) = Mijy
% s(4) = cos( theta_max )
% s(5) = theta_max
% s(6) = Dphi % orientation angle of vxi-vxj in global coords
% or, if wi=wj, a contact plane
% s(1) = NaN
% s(2) = x0  % starting point of
% s(3) = y0  % straight line
% s(4) = dx  % direction of 
% s(5) = dy  % straight line
% s(6) = Dphi
%
% each cell organized as
% c(1) = x
% c(2) = y
% c(3) = r
% c(4) = w
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s = [ s1; s2; s3 ];
sc = [ 2;  3;  1 ];    % cyclic


% find intersection points
ptm = 0;
pt = NaN*ones(6,2);
for m=1:3
  mm = sc(m);
  Ra = s(m, 1);
  Rb = s(mm,1);

  if isnan(Ra) && isnan(Rb)
    va = s(m, 2:3);
    vb = s(mm,2:3);
    ua = s(m, 4:5);
    ub = s(mm,4:5);
    pts = intlip( va,ua, vb,ub, mindist );
    sizpts = size(pts);
    if sizpts(1) == 2
      error('stop vertex.m: coinciding straight lines unexpected');
    end
  elseif isnan(Ra)
    va = s(m, 2:3);
    ua = s(m, 4:5);
    Mb = s(mm,2:3);
    pts = intlis( va,ua, Mb,Rb, mindist );
    sizpts = size(pts);
  elseif isnan(Rb)
    vb = s(mm,2:3);
    ub = s(mm,4:5);
    Ma = s(m, 2:3);
    pts = intlis( vb,ub, Ma,Ra, mindist );
    sizpts = size(pts);
  else
    Ma = s(m, 2:3);
    Mb = s(mm,2:3);
    pts = intsph( Ma,Ra, Mb,Rb, mindist );
    sizpts = size(pts);
    if sizpts(1) == 3 % spheres coinciding
      error('stop vertex.m: coinciding spheres unexpected');
    end
  end

  sizpts = size(pts);
  for n = 1:sizpts(1)
    ptm = ptm+1;
    pt(ptm,:) = pts(n,:);
  end
end


% find tripels of equal points
p3  = zeros(2,2);
p3m = 0;
for p=1:ptm
  nip = 1; % number of identical points
  ptp = pt(p,:);
  %fprintf( 1, 'vcand = [ %.16f, %.16f ]\n', ptp(1),ptp(2) );
  for pp=p+1:ptm
    ptpp = pt(pp,:);
    vdist = ptp - ptpp;
    dist = sqrt( vdist*vdist' );
    if dist < mindist
      nip = nip + 1;
    end
  end
  if nip==3
    p3m = p3m + 1;
    p3(p3m,:) = ptp;
  end
end


% choose point
if p3m==2
  retval = p3;
elseif p3m==1
  retval = p3(1,:);
else
  %warning('stop vertex.m: no vertex found');
  %retval = [ 0, 0 ];
  % error in full simulation
  %error('stop vertex.m: no vertex found');
  retval = [];
end

