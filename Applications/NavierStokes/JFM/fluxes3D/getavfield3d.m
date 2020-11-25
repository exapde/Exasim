function avField = getavfield3d(udg,odg,param,porder)
gam = param(1);
gam1 = gam - 1.0;
rmin = param(5);
avbulk = param(6);
Hmin = 1.e-4;
avs = param(8);
avb0 = param(12);

% slope regularization parameter
alpha = 1.0e6;
% regularization parameters for the bulk viscosity
kb = avbulk;      % kb = 1.5 for Ma<6
sb0   = avb0;
sbmax = 5.0 / sqrt(gam*gam - 1.0);
sbmin = 0.0;
% regularization parameters for the shear viscosity
km = avs;
sm0   = 0.01;
smmax = 0.1;
smmin = 0.0;

% Get Metric terms
avfix   = odg(5);
hm      = odg(6);

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

S_xx = 0.5 * (ux + ux);
S_xy = 0.5 * (uy + vx);
S_xz = 0.5 * (uz + wx);
S_yx = 0.5 * (vx + uy);
S_yy = 0.5 * (vy + vy);
S_yz = 0.5 * (vz + wy);
S_zx = 0.5 * (wx + uz);
S_zy = 0.5 * (wy + vz);
S_zz = 0.5 * (wz + wz);
shear = sqrt(S_xx*S_xx + S_yx*S_yx + S_zx*S_zx + ...
             S_xy*S_xy + S_yy*S_yy + S_zy*S_zy + ...
             S_xz*S_xz + S_yz*S_yz + S_zz*S_zz);

% limit 
sigm = 1.0e3;
div_v = limiting(div_v,-sigm,sigm,alpha,-sigm);
vort = limiting(vort,0.0,sigm,alpha,0.0);
shear = limiting(shear,0.0,sigm,alpha,0.0);

DucrosRatio = div_v.*div_v ./ (div_v.*div_v + vort.*vort + 1.0e-16);
%ShearRatio = div_v.*div_v ./ (div_v.*div_v + shear.*shear + 1.0e-16);

% Dilatation Sensor sb
sb = - (hm./porder) .* (div_v./c_star) .* DucrosRatio; %.*ShearRatio;
sb = limiting(sb,sbmin,sbmax,alpha,sb0);
% Artificial Bulk viscosity
avb = r .* (kb.*hm./(porder)) .* sqrt(uv.*uv + vv.*vv + c_star.*c_star) .* sb;

% shear viscosity
DucrosRatio = shear.*shear./(16*div_v.*div_v + shear.*shear + 1.0e-16);
sm = (hm./porder).*shear.*DucrosRatio./c_star;
sm = limiting(sm,smmin,smmax,alpha,sm0);
avs = r .* (km.*hm./(porder)) .* sqrt(uv.*uv + vv.*vv + c_star.*c_star) .* sm;

% Assign artificial viscosities
avField(1) = avfix.*avb;
avField(2) = 0;
avField(3) = avfix.*avb;
avField(4) = avfix.*avb;






