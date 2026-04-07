function [wall,land] = lbseawall()

% Recreates the 2D cross-section of the Low Battery seawall (Charleston)
% in the same style/shape as the provided figure.

%% ---------------- Units ----------------
ft2m = 0.3048;

%% ---------------- Elevations (feet, from figure) ----------------
z_mud  = 0.00;
z_top  = 8.20;

%% ---------------- Geometry knobs (tuned to match the figure look) ----------------
% These are "shape controls" chosen to reproduce the pictured outline.
% If you want a slightly different visual match, adjust these.

fac = 1.5;
% Wall thickness / step offsets (feet)
x_step1_in    = 1.7/fac;     % first inset at the top
x_step2_in    = 2.7/fac;     % second inset
x_mid_out     = 3.7/fac;     % outward bulge at mid-height
x_lower_in    = 5.4/fac;     % inset near the lower step
x_toe_land    = 7.0/fac;     % toe landward x at mudline (where wall meets mud)

% Heights for intermediate steps (feet) (picked to match the drawing)
z1 = 6.50;              % small drop under cap
z2 = 4.80;              % sidewalk elevation (exact)
z3 = 3.30;              % mid step (schematic)
z4 = 1.80;              % lower step (schematic)

% Seaward face curve control (feet): (x,z) points (increasing z)
sea_curve = [
   -2.6  -1.0
   -1.8  -0.5
   -1.2  -0.2
   -0.7   1.2
   -0.4   3.2
   -0.2   5.3
   -0.1   6.8
    0.0   z_top
];

% Stepped profile (matches the figure’s “stair-step” interior)
land = [
    x_step1_in  z_top
    x_step1_in  z1
    x_step2_in  z1
    x_step2_in  z2
    x_mid_out   z2
    x_mid_out   z3
    x_lower_in  z3
    x_lower_in  z4
    x_toe_land  z4
    x_toe_land  z_mud    
];

% Seaward face: spline-interpolate the given curve control points
zq = linspace(sea_curve(1,2), sea_curve(end,2), 250).';
xq = interp1(sea_curve(:,2), sea_curve(:,1), zq, 'pchip');  % x(z)

sea = [xq, zq];  % from bottom up to crest on seaward side
see = sea(sea(:,2)>=z_mud,:);
see = [-1.074154005128054  z_mud; see]; 
a = [see(1:32:end,:); see(end,:)];
a(end-1,:) = [];

wall = [a; land]*ft2m;
land = land*ft2m;
