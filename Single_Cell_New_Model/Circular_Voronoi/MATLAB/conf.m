%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Plotting setup
asprat = 4/3;          % aspect ratio for plotting

plot_cell_centnu = 1;  % plot cell center  numbers        1:y else no
plot_cell_body   = 1;  % plot cell bodies                 1:y else no
plot_contact_arc = 1;  % plot contact arcs                1:y else no
plot_margin_arc  = 1;  % plot marginal arcs               1:y else no
plot_delaunay    = 1;  % plot Delaunay triang             1:y else no
plot_axis        = 1;  % plot axis                        1:y else no

col_cell = 'g';        % color of cell                    green
col_contact_arc = 'r'; % color of contact arcs            red
col_margin_arc  = 'k'; % color of marginal arcs           black
col_delaunay = 'g';    % color of delaunay triangulation  green

lw_cell_body   = 1;    % linewidth cell body
lw_contact_arc = 1;    % linewidth contact arc

npts_contact_arc = 064; % number of points to approximate cell-cell contacts
npts_margin_arc  = 064; % number of points to approximate cell closure circles
npts_cell_body   = 032; % number of points to approximate cell body circles

static_scale = 0;    % apply static plot scaling
statscal = [ -40,40, -30,30 ];





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% program flow
runfromoct = 0;      % run from octave?

% what kind of configurations shall be allowed
check_thmax = 0;     % check wether cells are starlike (theta < thmax)
check_digestion = 1; % check wether cells are digesting 
                     % (overlap w/o pseudo wertices)

% Pmax, and equality thresholds
Pmax = 3.00;           % maximal distance function sqrt(P_max)
mindist = 1e-6;        % minimal distance, below points are considered equal
minang  = 1e-6;        % minimal distance, below angles are considered equal

% allocation constants
nol = 32;              % max # of possible overlaps to assume for single cell
                       % utmost upper bound, probably much too high
                       % also used as estimate for various other index maxima

% paths
addpath('./func');





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% cell data
typ = [ 1,  1,     1,    1,     1,     1,     1,   ]; % type of cell
                                                      % no significance here

x = [ 1.00, 1.00,  3.00, 6.00, 10.00,  2.50, -5.00 ]; % x positions nuclei
y = [ 4.00, 0.00, 10.00, 0.00,  9.00, -5.00,  8.00 ]; % y positions nuclei
w = [ 1.50, 1.30,  1.85, 1.70,  1.75,  1.40,  1.60 ]; % weights
r = w;

cm = length(x);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% building tessellation
%toggle_octplot;
fig = figure();


mwvoro;
plotcells;
input('press enter to continue');

