%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% mwvoro.m
%%%%%%%% compute circular voronoi tesselation
%%%%%%%%
%%%%%%%% License:
% Copyright (c) 2009
% Martin Bock. All rights reserved.
%
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
%    this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation 
%    and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE AUTHOR "AS IS" AND ANY EXPRESS OR IMPLIED 
% WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO 
% EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
% SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
% OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
% WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
% OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
% ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% short explanation of the spirit of variable naming:
%
% c Cell
% p (cell) Pair
% t (cell) Tripel (pairwise overlapping)
% a contact arc candidate
% f free closure arc
%
% i,j,k,l  % cell indices
% m,mm,... % overlapping cell pair indices
% n,nn,... % overlap triple (or vertex) indices
% o,oo,... % pseudo vertex indices
% p,pp,... % Voronoi arc indices
% q,qq,... % free, marginal arc indices


% relating and linking the various involved objects
% Cell -> Pairs, Triples
%cm                        % # cells
c2a = zeros(cm,nol);       % indices of Voronoi arcs arround cell {p}
c2p = zeros(cm,nol);       % indices of pairs containing cell {m}
c2am = zeros(1,cm);        % # of Voronoi arcs arround cell
c2pm = zeros(1,cm);        % # of pairs containing cell


% Pair -> Cells, Triples, Vertices
pm = 0;                    % # pairs
p2s = zeros(nol*cm,7);     % contact sphere  1:Rij  2:Mijx  3:Mijy
                           %   4:cos(thmax)  5:thmax
                           %   6:thij  orientation angle of vx_big - vx_small
                           %   7:delta cell body surface distance
                           % or plane: 1:NaN  2:x0   3:y0  (local origin)
                           %   4:slpx 5:slpy  normalized slope of straight line
                           %   6: thij  orientation angle of vx_i - vx_j
                           %   7:delta cell body surface distance
p2c = zeros(nol*cm,2);     % indices of cells i,j forming pair
p2t = zeros(nol*cm,nol);   % indices of tripels containing pair n,nn,...
p2v = zeros(nol*cm,nol);   % vertices of pair n,nn,...
p2tm = zeros(1,nol*cm);    % # tripels containing pair
p2vm = zeros(1,nol*cm);    % # vertices of pair


% Tripel -> Cells, Pairs
tm = 0;                    % # triples
t2c = zeros(nol*cm,3);     % indices of cells forming tripel i,j,k
t2p = zeros(nol*cm,3);     % indices of pairs forming tripel m,mm,...


% Vertices, Vertices -> Cells
vm = 0;                    % # vertices
v2v = zeros(nol*cm,3);     % vertex coords 1:x 2:y, 3:P_distance
v2c = zeros(nol*cm,3);     % indices of cells forming vertex i,j,k


% pseudo Wertices arising from Pmax
wm = 0;                    % # pseudo vertices
w2w = zeros(2*nol*cm,3);   % pseudo vertex coordinates, 1:x 2:y, 3:P_distance
p2w = zeros(nol*cm,2);     % pseudo vertices of pair o,oo


% contact arcs
am = 0;                    % number of contact arc candidates
a2a = zeros(nol*cm,14);    % 1:m  pair index
                           % 2:R  radius of contact sphere
                           % 3:Mx 4:My  center of sphere
                           % 5:theta1 6:theta2  angle range of contact surface
                           % 7,8 + 9,10:     start and end vertex
                           % 11,12  type  of start and and vertex 0w/1t
                           % 13,14  index of start and and vertex  o/n
                           % or flat contact surface
                           % 2: NaN
                           % 3:xi 4:yi  coords of cell i
                           % 5:0  6:0


% free arcs
fm = 0;                    % number of free closure arcs
f2f = zeros(nol*cm,6);     % 1: cell index i  2: Pmax*w(i)
                           % 3:xi 4:yi  cell center
                           % 5:thmin 6:thmax w.r.t. cell center


% Voronoi/Delaunay neighbors
nb   = zeros(cm,nol);      % indices of Voronoi neighbors of cell
nbm  = zeros(1,cm);        % number of Voronoi neighbors of cell
nbd  = zeros(1,nol*cm);    % neighbor pair
nbdm = 0;                  % # neighbor pairs


% free cell array
fcm = 0;                   % # of free cells
fc  = zeros(cm,1);         % indices of free cells


% marginal cell array
mcm = 0;                   % # of marginal cells
mc  = zeros(cm,1);         % indices of marginal cells


% % temporary variables for voronoi partition
vcl  = zeros(2*nol,6);     % vertex candidate list per pair
                           % 1:orientation angle thv towards Mij or cell i
                           % 1:P_distance w.r.t. cells,
                           % 2:x 3:y  vertex coordinates
                           % 4:index n,o  5:type=0w/1t/2u of v/w/uertex
vl   = zeros(2*nol,6);     % vertex list per pair
                           % 1:orientation angle thv (towards Mij or left cell),
                           % 2:P_distance w.r.t. first cell,
                           % 3:x 4:y  vertex coordinates
                           % 5:index n,o  6:type=0w/1t/2u of v/wertex
svl  = zeros(2*nol,6);     % sorted vertex list per pair
                           % 1:orientation angle thv (towards Mij or left cell),
                           % 2:P_distance w.r.t. first cell,
                           % 3:x 4:y  vertex coordinates
                           % 5:index n,o  6:type=0w/1t/2u of v/wertex
al  = zeros(nol,14);       % per cell arc list, cf. a2a above, 13:phi1 14:phi2
sal = zeros(nol,14);       % per cell arc list, sorted after phi1, the
                           % angular starting point of the arc w.r.t. cell
                           % and phi2 angular ending point of arc w.r.t. cell





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% determine pairs of overlapping cells
for i = 1:cm
  jj = i+1:cm;
  distx = x(i) - x(jj);
  disty = y(i) - y(jj);
  dist = sqrt( distx.*distx + disty.*disty );
  cclose = i + find(  dist < Pmax*( w(i) + w(jj) )  );
  sizcl = size(cclose);
  if sizcl(1) > 0
    for j = cclose
      pm = pm + 1;
      m  = pm;
      % i,j sorted in ascending order
      p2c(m,:) = [ i, j ];
      c2pm(i) = c2pm(i) + 1;
      c2pm(j) = c2pm(j) + 1;
      c2p( i,   c2pm(i) ) = m;
      c2p( j,   c2pm(j) ) = m;
    end
  end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% calculate contact spheres
for m = 1:pm % loop over pairs

  i = p2c(m,1);
  j = p2c(m,2);
  if abs( w(i)-w(j) ) < mindist
    sgn =    1;
    cthm = w(j)/w(i);
  elseif w(i) > w(j)
    sgn  =   1;
    cthm = w(j)/w(i);
  else
    sgn  = - 1;
    cthm = w(i)/w(j);
  end
  dx = ( x(i) - x(j) )*sgn;
  dy = ( y(i) - y(j) )*sgn;
  x0 = ( x(i) + x(j) )/2;
  y0 = ( y(i) + y(j) )/2;
  dij = sqrt( dx^2 + dy^2 );
  delta = dij - r(i) - r(j);

  thij = atan2( dy, dx );
  if abs( w(i)-w(j) ) < mindist
    % direction of straight line starting from pair-local origin 
    % (x0,y0) given by (-dy,dx) in global coords
    vdist = [ dx, dy ];  % record orientation of pair
    udist = vdist/sqrt( vdist*vdist' );
    slp = [ - udist(2), udist(1) ];
    Rij = NaN;
    p2s(m,:) = [ Rij, x0,y0, slp(1),slp(2), thij, delta ];
  else
    % the usual contact sphere
    thm = acos(cthm);        % numerically unstable at wi ~ wj
    wwp =     w(i)^2 + w(j)^2;
    wwm = ( - w(i)^2 + w(j)^2 )*sgn;
    Rij = dij*w(i)*w(j)*sgn/( w(i)^2-w(j)^2 );
    vMij = ShiftRot2d( wwp*dij/wwm/2,0, x0,y0, thij );
    Mijx = vMij(1);
    Mijy = vMij(2);
    p2s(m,:) = [ Rij, Mijx,Mijy, cthm,thm, thij, delta ];
  end

end % loop over pairs





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% determine tripels of overlapping cells and vertices 
for m = 1:pm % loop over pairs m
  for mm = m+1:pm % loop over pairs mm
    % due to sortedness comparing these indices is sufficient
    if p2c(m,2) == p2c(mm,1) % if pairs share common cell
      i = p2c( m,  1 ); % i,j,k sorted in ascending order
      j = p2c( mm, 1 ); % pair  m = i,j
      k = p2c( mm, 2 ); % pair mm = j,k
      m3 = intersect( c2p(i,1:c2pm(i)), c2p(k,1:c2pm(k)) ); % pair m3 = i,k
      sizm3 = size(m3); % sizm3 = [1 1] or [0 0], could be checked
      if sizm3 == [1,1] % check if actual tripel
        ii = p2c(m3,1);
        kk = p2c(m3,2);
        if ii ~= i || kk ~= k
          warning('warn inivor.m: strange tripel rejected');
          continue;
        end
      else
        continue;
      end

      % associating overlap tripel
      tm = tm + 1;
      n  = tm;
      t2p(n,:) = [ m, m3, mm ]; % indeed sorted due to loop construction
      t2c(n,:) = [ i, j, k ];
      p2tm( m) = p2tm( m) + 1;
      p2tm(mm) = p2tm(mm) + 1;
      p2tm(m3) = p2tm(m3) + 1;
      p2t(  m, p2tm( m) ) = n;
      p2t( mm, p2tm(mm) ) = n;
      p2t( m3, p2tm(m3) ) = n;

      % associating vertex
      ci = [ x(i), y(i), r(i), w(i) ];
      cj = [ x(j), y(j), r(j), w(j) ];
      ck = [ x(k), y(k), r(k), w(k) ];
      vrtx = vertex( p2s(m,:), p2s(m3,:), p2s(mm,:), ci,cj,ck, mindist );
      sizvrtx = size(vrtx); % 0, 1 or 2 points
      for nn = 1:sizvrtx(1)
        vdisti = [x(i),y(i)]-vrtx(nn,:);
        vdistj = [x(j),y(j)]-vrtx(nn,:);
        vdistk = [x(k),y(k)]-vrtx(nn,:);
        Pi = sqrt( vdisti*vdisti' )/w(i);
        Pj = sqrt( vdistj*vdistj' )/w(j);
        Pk = sqrt( vdistk*vdistk' )/w(k);
        if Pi < Pmax  &&  Pj < Pmax  && Pk < Pmax
          vm = vm + 1;
          v2v(vm,:) = [ vrtx(nn,1),vrtx(nn,2), Pi ];
          v2c(vm,:) = [ i, j, k ];
          p2vm( m) = p2vm( m) + 1;
          p2vm(mm) = p2vm(mm) + 1;
          p2vm(m3) = p2vm(m3) + 1;
          p2v(  m, p2vm( m) ) = vm;
          p2v( mm, p2vm(mm) ) = vm;
          p2v( m3, p2vm(m3) ) = vm;
        end
      end

    end % if pairs share common cell
  end % loop over pairs mm
end % loop over pairs m



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% calculate pseudo-vertices of cell pairs
%%% Due to limiting the extension of cells up to the maximal 
%%% generalized distance Pmax, for each overlapping cell pair there
%%% appear two additional pseudo-vertices given by the intersection 
%%% of both cell's Pmax-spheres.
for m = 1:pm % loop over pairs
  i = p2c(m,1);
  j = p2c(m,2);
  si = [ Pmax*w(i), x(i),y(i) ];
  sj = [ Pmax*w(j), x(j),y(j) ];
  pts = wertex( si, sj, mindist );
  sizpts = size(pts);
  if check_digestion == 1 && sizpts(2) < 2
    error('stop inivor.m: cells overlapping w/o wertex, digestion detected');
  end
  if sizpts(1) == 1 % singular point overlap ignored
    continue;
  end
  for o=1:sizpts(1)
    pto = pts(o,:);
    wm = wm + 1;
    w2w(wm,:) = [ pto(1),pto(2), Pmax ];
    p2w( m,o) = wm;
  end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% construct interior Voronoi arcs and neighbor relations
for m = 1:pm % loop over pairs

  % initial variables, zero counters
  i = p2c(m,1);
  j = p2c(m,2);
  %printf('pair %i, cells %i,%i\n', m, i, j );
  Rij = p2s(m,1);
  thmax = p2s(m,5);
  thij  = p2s(m,6);
  vclm  = 0; % # vertex candidates of pair
  vlm   = 0; % # w/vertices of pair
  neighbor = 0; % is pair neighbor pair?


  % build neighbor candidate list of pair i,j from their tripels
  nn = p2t(m,1:p2tm(m));
  kk = vec(t2c(nn,:))';
  kk = unique(kk);
  if runfromoct == 1
    cnbc = complement( [i,j], kk );    % cnbc common neighbor candidates
  else
    cnbc = setdiff( kk, [i,j] );       % order indeed reversed
  end


  % build vertex candidate list (trading redundancy in code for performance)
  if isnan(Rij) % flat contact surface
    vx0 = p2s(m,2:3);
    for no=1:p2vm(m) % vertices
      vclm = vclm + 1;
      n = p2v(m,no);
      vrtx  = v2v(n,1:2);
      Pvrtx = v2v(n, 3 );
      vcl(vclm,2:6) = [ Pvrtx, vrtx(1),vrtx(2), n,1 ];
    end
    for no=1:2 % pseudo wertices
      o = p2w(m,no);
      if o ~= 0 % checking for digestion earlier
        vclm = vclm + 1;
        wrtx  = w2w(o,1:2);
        Pwrtx = w2w(o, 3 ); % = Pmax
        vcl(vclm,2:6) = [ Pwrtx, wrtx(1),wrtx(2), o,0 ];
      end
    end
    % The artificial u-vertex is needed to decide
    % wether the arc is to the left or right of Mij.
    % By separating at the critical angle theta = +/- pi 
    % (think atan2), all angle ranges are contiguous.
    % Later on this thingie shall be deleted, unless it forms
    % the only vertex pair of the arc, like in digesting systems.
    if abs( x(i) - x(j) ) > mindist && x(i) > x(j)
      ux0 = p2s(m,4:5);
      urtx = intlip( vx0,ux0, [x(i),y(i)],[1,0], mindist);
      sizurtx = size(urtx);
      if sizurtx(1) == 1
        vdist = [x(i),y(i)] - urtx;
        Purtx = sqrt( vdist*vdist' )/w(i);
        if Purtx < Pmax
          vclm = vclm + 1;
          vcl(vclm,:) = [ -pi, Purtx, urtx(1),urtx(2), no,2 ];
          vclm = vclm + 1;
          vcl(vclm,:) = [  pi, Purtx, urtx(1),urtx(2), no,2 ];
        end
      end
    end
  else % spherical contact surface
    vMij = p2s(m,2:3);
    for no=1:p2vm(m)
      vclm = vclm + 1;
      n = p2v(m,no);
      vrtx  = v2v(n,1:2);
      Pvrtx = v2v(n, 3 );
      vcl(vclm,2:6) = [ Pvrtx, vrtx(1),vrtx(2), n,1 ];
    end
    for no=1:2
      o = p2w(m,no);
      if o ~= 0 % checking for digestion earlier
        vclm = vclm + 1;
        wrtx  = w2w(o,1:2);
        Pwrtx = w2w(o, 3 ); % = Pmax
        vcl(vclm,2:6) = [ Pwrtx, wrtx(1),wrtx(2), o,0 ];
      end
    end
    % cf. urtx comment at planar contact above
    urtx = vMij - Rij*[1,0];
    vdist = [x(i),y(i)] - urtx;
    Purtx = sqrt( vdist*vdist' )/w(i);
    if Purtx < Pmax
      vclm = vclm + 1;
      vcl(vclm,:) = [ -pi, Purtx, urtx(1),urtx(2), no,2 ];
      vclm = vclm + 1;
      vcl(vclm,:) = [  pi, Purtx, urtx(1),urtx(2), no,2 ];
    end
  end % if flat contact surface


  % vertex erosion: find vertices which are not closer to other cells
  for no = 1:vclm
    Pi   = vcl(no,2);
    vrtx = vcl(no,3:4);
    nnoo = vcl(no,5);
    vtyp = vcl(no,6);
    if vtyp == 1 % avoid comparison of equal P_distances
      ijk = v2c(nnoo,:);
      if runfromoct == 1
        ccnbc = complement( ijk, cnbc ); % cnbc common neighbor candidates
      else
        ccnbc = setdiff( cnbc, ijk ); % order indeed reversed
      end
    else
      ccnbc = cnbc;
    end
    bbreak = 0;
    for k = ccnbc
      vdistk = [x(k),y(k)]-vrtx;
      Pk = sqrt(vdistk*vdistk')/w(k);
      if Pk < Pi % vertex closer to other cell
        bbreak = 1;
        break;
      end
    end
    if bbreak == 1 % disregard vertex
      continue;
    else % associate vertex
      if isnan(Rij) 
        if vtyp == 0  ||  vtyp == 1
          vdist = vrtx - [x(i),y(i)];
          vcl(no,1) = atan2( vdist(2), vdist(1) );
        end
      else
        if vtyp == 0  ||  vtyp == 1
          vdist = vrtx - vMij;
          vcl(no,1) = atan2( vdist(2), vdist(1) );
        end
      end
      vlm = vlm + 1;
      vl(vlm,:) = vcl(no,:);
    end
  end
  

  % remove multiply appearing vertices for degenerate case
  if mod(vlm,2) ~= 0
    vclm = 0;
    for no = 1:vlm
      vno = vl(no,3:4);
      vtyp = vl(no,6);
      if vtyp == 2 % urtx always needed
        vclm = vclm + 1;
        vcl(vclm,:) = vl(no,:);
      else
        solo = 1;
        for on = no+1:vlm
          von = vl(on,3:4);
          vdv = vno - von;
          dv = sqrt( vdv*vdv' );
          if dv < mindist
            solo = 0;
          end
        end
        if solo == 1;
          vclm = vclm + 1;
          vcl(vclm,:) = vl(no,:);
        end
      end
    end
    vlm = vclm;
    vl(1:vlm,:) = vcl(1:vclm,:);
  end
  [vvl,no2] = sort( vl(1:vlm,1) ); % sort by thv
  svl(1:vlm,:) = vl(no2,:);


  % determine and associate arcs
  if mod(vlm,2) ~= 0
    error('stop inivor.m: odd number of vertices at i=%i, j=%i',i,j);
    %warning('stop inivor.m: odd number of vertices at i=%i, j=%i',i,j);
    %input('press enter to continue');
  end
  vM = p2s(m,2:3);
  vtyp = svl(1,6);
  if vlm >= 4 && vtyp == 2 % first reconnect splitted arc
    th1 = svl( vlm-1,  1  ); % since first and last uertex are pseudo
    th2 = svl(   2  ,  1  ) + 2*pi;
    v1   = svl( vlm-1, 3:4 );
    v2   = svl(   2,   3:4 );
    v1no = svl( vlm-1,  5  );
    v2no = svl(   2,    5  );
    v1t  = svl( vlm-1,  6  );
    v2t  = svl(   2,    6  );
    if check_thmax == 1
      dth1 = th1 - thij;
      dth2 = th2 - thij;
      if abs(dth1) > thmax || abs(dth2) > thmax
        error('stop inivor.m: starlikeness violated');
      end
    end
    neighbor = 1;
    am = am + 1;
    a2a(am,:) = [ m, Rij, vM(1),vM(2), th1,th2, v1(1),v1(2), v2(1),v2(2), ...
                                                  v1t,v2t,    v1no,v2no       ];
    c2am(i) = c2am(i) + 1;
    c2am(j) = c2am(j) + 1;
    c2a( i,   c2am(i) ) = am;
    c2a( j,   c2am(j) ) = am;
    nbm(i) = nbm(i) + 1;
    nbm(j) = nbm(j) + 1;
    nb( i,   nbm(i) ) = j;
    nb( j,   nbm(j) ) = i;
    noran = 3:2:vlm-2; % so uertex is deleted
  else
    noran = 1:2:vlm;
  end
  for no=noran % then proceed to associate usual arcs
    v1   = svl(no  ,3:4);
    v2   = svl(no+1,3:4);
    v1no = svl(no,   5 );
    v2no = svl(no+1, 5 );
    v1t  = svl(no  , 6 );
    v2t  = svl(no+1, 6 );
    if isnan(Rij)
      th1 = thij;
      th2 = thij;
    else
      th1 = svl(no  ,1);
      th2 = svl(no+1,1);
    end
    if check_thmax == 1
      dth1 = th1 - thij;
      dth2 = th2 - thij;
      if abs(dth1) > thmax || abs(dth2) > thmax
        error('stop inivor.m: starlikeness violated');
      end
    end
    neighbor = 1;
    am = am + 1;
    % Rij = Nan if planar contact
    a2a(am,:) = [ m, Rij, vM(1),vM(2), th1,th2, v1(1),v1(2), v2(1),v2(2), ...
                                                  v1t,v2t,    v1no,v2no       ];
    c2am(i) = c2am(i) + 1;
    c2am(j) = c2am(j) + 1;
    c2a( i,   c2am(i) ) = am;
    c2a( j,   c2am(j) ) = am;
    nbm(i) = nbm(i) + 1;
    nbm(j) = nbm(j) + 1;
    nb( i,   nbm(i) ) = j;
    nb( j,   nbm(j) ) = i;
  end
  if neighbor == 1
    nbdm = nbdm + 1;
    nbd(nbdm) = m;
  end

end % loop over pairs



% find marginal Voronoi arcs comprising the free boundary
for i = 1:cm
  atmargin = 0;
  alm = 0;
  numdig = 0; % number of digested cells
  posi = [x(i),y(i)]; % cell center
  for p = 1:c2am(i)
    alc = a2a( c2a(i,p), : ); % arc list candidate
    m = alc(1);
    ij = p2c(m,:);
    if runfromoct == 1
      j  = complement( [i], ij );
    else
      j  = setdiff( ij, [i] ); % order indeed reversed
    end
    Ri = Pmax*w(i);
    Rj = Pmax*w(j);
    posj = [x(j),y(j)];
    vdij = posi - posj;
    dij  = sqrt( vdij*vdij' );
    if Ri-Rj > dij % cell i digesting cell j 
      numdig = numdig + 1;
    end
    if alc(11) == 0 || alc(12) == 0 % consider only arcs w/ wertices
      if abs( w(i) - w(j) ) < mindist
        if i > j % vertices have been sorted w.r.t. first cell of pair
          v1  = alc(9:10);
          v2  = alc(7:8);
          vt1 = alc(12);
          vt2 = alc(11);
        else
          v1  = alc(7:8);
          v2  = alc(9:10);
          vt1 = alc(11);
          vt2 = alc(12);
        end
      else
        if w(i) > w(j)
          v1  = alc(9:10);
          v2  = alc(7:8);
          vt1 = alc(12);
          vt2 = alc(11);
        else
          v1  = alc(7:8);
          v2  = alc(9:10);
          vt1 = alc(11);
          vt2 = alc(12);
        end
      end
      vdist1 = v1 - posi;
      vdist2 = v2 - posi;
      phi1 = atan2( vdist1(2), vdist1(1) );  % orientation starting vertex
      phi2 = atan2( vdist2(2), vdist2(1) );  % orientation  ending  vertex
      alm = alm + 1;
      al(alm,:) = [ alc(1:6), v1,v2, vt1,vt2, phi1,phi2 ];
    end
  end
  [aal,pp] = sort( al(1:alm,13) );
  sal(1:alm,:) = al(pp,:);
  if alm == 0
    if numdig - c2am(i) == 0 % free or purely digesting cell, complete circle
      fm = fm + 1;
      f2f(fm,:) = [ i, Pmax*w(i), x(i),y(i), -pi,pi ];
      fcm = fcm + 1;
      fc(fcm) = i;
      atmargin = 1;
    end
  else % connect wertices for cells w/ contact arcs
    for p = 1:alm-1
      pp = p + 1;
      vt1 = sal(p ,12);
      vt2 = sal(pp,11);
      if vt1==0 && vt2==0
        phi1 = sal(p ,14);
        phi2 = sal(pp,13);
        fm = fm + 1;
        f2f(fm,:) = [ i, Pmax*w(i), x(i),y(i), phi1,phi2 ];
        atmargin = 1;
      end
    end
    vt1 = sal(alm,12);
    vt2 = sal( 1 ,11);
    if vt1==0 && vt2==0
      phi1 = sal(alm,14);
      phi2 = sal( 1 ,13);
      if phi2 < phi1
        phi2 = phi2 + 2*pi;
      end
      fm = fm + 1;
      f2f(fm,:) = [ i, Pmax*w(i), x(i),y(i), phi1,phi2 ];
      atmargin = 1;
    end
  end
  if atmargin == 1
    mcm = mcm + 1;
    mc(mcm) = i;
  end
end

