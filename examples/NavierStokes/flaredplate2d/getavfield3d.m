function avField = getavfield3d(udg,odg,param,porder)

% odg = v
% v(1) = AV field
% v(2) = mesh size field
% v(3) = laminar density field
% v(4:6) = laminar momentum field
% v(7) = laminar energy field 


% v(1) = laminar density field
% v(2:4) = laminar momentum field
% v(5) = laminar energy field 
% v(6) = AV field

% v(1) -> AV field 

% mesh size
hm = odg(2);

gam = param(1);
gam1 = gam - 1.0;

% artificial viscosity
avbulk = param(7);
avc = param(10);

% regularization parameters
alpha = 1.0e3;
rmin = 1.0e-3;
Hmin = 1.0e-4;

% regularization parameters for the bulk viscosity
kb = avbulk;      % kb = 1.5 for Ma<6
sb0   = 0.05;
sbmax = 5.0 / sqrt(gam*gam - 1.0);
sbmin = 0.0;

% Get base variables
r = udg(1);
ru = udg(2);
rv = udg(3);
rw = udg(4);
rE = udg(5);
rx = udg(6);
rux = udg(7);
rvx = udg(8);
rwx = udg(9);
%rEx = udg(10);
ry = udg(11);
ruy = udg(12);
rvy = udg(13);
rwy = udg(14);
%rEy = udg(15);
rz = udg(16);
ruz = udg(17);
rvz = udg(18);
rwz = udg(19);
%rEz = udg(20);

% Regularization of density 
r = rmin + lmax(r-rmin,alpha);
r1 = 1./r;
uv = ru.*r1;
vv = rv.*r1;
wv = rw*r1;
E = rE.*r1;
q = 0.5*(uv*uv+vv*vv+wv*wv);
H = gam*E - gam1*q;                    % Enthalpy (Critical Speed of Sound ???)
H = Hmin + lmax(H-Hmin,alpha);         % Regularized Enthalpy

% Critical speed of Sound
c_star = sqrt((2.*gam1.*H) ./ (gam+1));

% Computing derivatives for the sensors
ux = (rux - rx*uv)*r1;
vx = (rvx - rx*vv)*r1;
wx = (rwx - rx*wv)*r1;
uy = (ruy - ry*uv)*r1;
vy = (rvy - ry*vv)*r1;
wy = (rwy - ry*wv)*r1;
uz = (ruz - rz*uv)*r1;
vz = (rvz - rz*vv)*r1;
wz = (rwz - rz*wv)*r1;

div_v = - (ux + vy + wz);
vort_x = - (wy - vz);
vort_y = - (uz - wx);
vort_z = - (vx - uy);
vort = sqrt(vort_x*vort_x + vort_y*vort_y + vort_z*vort_z);


% limit  divergence and vorticity
sigm = 1e3;
div_v = limiting(div_v,-sigm,sigm,alpha,-sigm);
vort = limiting(vort,0.0,sigm,alpha,0.0);

% Dilatation Sensor sb
DucrosRatio = div_v.*div_v ./ (div_v.*div_v + vort.*vort + 1.0e-16);
% DucrosRatio = 1.0;
sb = - (hm./porder) .* (div_v./c_star) .* DucrosRatio;
sb = limiting(sb,sbmin,sbmax,alpha,sb0);
% Artificial Bulk viscosity
avb = r.*(kb.*hm./(porder)) .* sqrt(uv.*uv + vv.*vv + c_star.*c_star) .* sb;

avb = avb;

% Assign artificial viscosities
avField(1) = avb;  %  bulk

