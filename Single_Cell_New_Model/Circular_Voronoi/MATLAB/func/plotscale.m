function retval = plotscale(xmin,xmax,ymin,ymax,radmax,asprat,fig)
% xmin: minimal x coord
% xmax: maximal x coord
% ymin: minimal y coord
% ymax: maximal y coord
% radmax: maximal radius arround coord point
% asprat: aspect ratio
% fig: figure

addmargin = 0.1; % how much margin to add in terms of radmax

rradm = (1+addmargin)*radmax;
xxmin = xmin - rradm;
xxmax = xmax + rradm;
yymin = ymin - rradm;
yymax = ymax + rradm;

ran(1)  = xxmax - xxmin;
ran(2)  = yymax - yymin;
yorx    = ( ran(2) > ran(1) );
xxmult  = ( 1 + (asprat-1)*(  yorx) )^( 2*(asprat>=1) - 1 );
yymult  = ( 1 + (asprat-1)*(1-yorx) )^( 2*(asprat>=1) - 1 );
ran(1)  = ran(1) / xxmult;
ran(2)  = ran(2) * yymult;
[maxxy coord] = max([ ran(1), ran(2) ]);
ran(1)  = ran(1) * xxmult;
ran(2)  = ran(2) / yymult;

ran(3-coord) = asprat^(2*coord-3)*ran(coord);

xav = xxmax+xxmin;
yav = yymax+yymin;
xxmin = ( xav - ran(1) )/2.0;
xxmax = ( xav + ran(1) )/2.0;
yymin = ( yav - ran(2) )/2.0;
yymax = ( yav + ran(2) )/2.0;

retval = [xxmin, xxmax, yymin, yymax];

