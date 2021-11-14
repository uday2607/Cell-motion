%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotcells.m
%
% plots cell configuration
% ATTENTION: all quantities defined and explained in main.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% scaling and preliminaries
clf;
if plot_axis == 0
  axis('off');
end
if runfromoct == 1
  rmax = max(r)*Pmax;
  if static_scale ~= 1
    xmin = min(x);
    xmax = max(x);
    ymin = min(y);
    ymax = max(y);
    rmax = max(r)*Pmax;
    winsiz = plotscale( xmin, xmax, ymin, ymax, rmax, asprat, fig );
  else
    axis( statscal );
    winsiz = statscal;
  end
else
  if static_scale ~= 1
    axis('equal');
  else
    winsiz = statscal;
  end
end
hold on;



% delaunay triangulation: linking of neighbor cells
if plot_delaunay == 1
  for mm = 1:nbdm
    m = nbd(mm);
    i = p2c(m,1);
    j = p2c(m,2);
    plx = [x(i),x(j)];
    ply = [y(i),y(j)];
    line( plx,ply, 'Color',col_delaunay );
  end
end



% cell bodies, Pmax spheres
if plot_cell_body == 1
  for i = 1:cm
    pl1 = arc( r(i), [x(i),y(i)], [-pi,pi], npts_cell_body, minang );
    line( pl1(1,:), pl1(2,:), 'Color',char( col_cell(typ(i)) ), ...
                              'LineWidth', lw_cell_body            );
  end
end



% cell-cell contacts
if plot_contact_arc == 1
  for p = 1:am
    Rij = a2a(p,2);
    if isnan( Rij )
      v1 = a2a(p,7:8);
      v2 = a2a(p,9:10);
      plx = [ v1(1), v2(1) ];
      ply = [ v1(2), v2(2) ];
      line( plx,ply, 'Color',col_contact_arc, 'LineWidth', lw_contact_arc );
    else
      Mij = a2a(p,3:4);
      thmin = a2a(p,5);
      thmax = a2a(p,6);
      pl = arc( Rij, Mij, [thmin,thmax], npts_contact_arc, minang );
      line( pl(1,:), pl(2,:), 'Color',col_contact_arc, ...
                              'LineWidth', lw_contact_arc  );
    end
  end
end



% free arcs of marginal cells
if plot_margin_arc == 1
  for q = 1:fm
    R   = f2f(q,2);
    vM  = f2f(q,3:4);
    phi = f2f(q,5:6);
    pl = arc( R, vM, phi, npts_margin_arc, minang );
    line( pl(1,:), pl(2,:), 'Color',col_margin_arc );
  end
end



% cell numbers
if plot_cell_centnu == 1
  for i = 1:cm
    cnum = num2str(i);
    text( x(i), y(i), cnum, 'Color', char( col_cell(typ(i)) ) );
  end
end



% finishing
if runfromoct==1
  axis(winsiz);
  octplot_command('redraw');
else
  drawnow;
end

