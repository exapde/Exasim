function f = getfluxav3d(udg,odg,param)
gam = param(1);
gam1 = gam - 1.0;
Re = param(2);
Pr = param(3);
Minf = param(4);
rmin = param(5);
pmin = param(7);
Tref = param(9);
alpha = 1.0e3;
muRef = 1/Re;
Tinf = 1/(gam*gam1*Minf^2);
c23 = 2.0/3.0;
avb = odg(1); % Bulk    viscosity
avr = odg(2); % Density viscosity
avs = odg(3); % Shear   viscosity
avk = odg(4); % Thermal viscosity

r = udg(1);
ru = udg(2);
rv = udg(3);
rw = udg(4);
rE = udg(5);
rx = udg(6);
rux = udg(7);
rvx = udg(8);
rwx = udg(9);
rEx = udg(10);
ry = udg(11);
ruy = udg(12);
rvy = udg(13);
rwy = udg(14);
rEy = udg(15);
rz = udg(16);
ruz = udg(17);
rvz = udg(18);
rwz = udg(19);
rEz = udg(20);

r = rmin + lmax(r-rmin,alpha);
dr = atan(alpha*(r - rmin))/pi + (alpha*(r - rmin))/(pi*(alpha^2*(r - rmin)^2 + 1)) + 1/2;
rx = rx*dr;
ry = ry*dr;
rz = rz*dr;

r1 = 1/r;
uv = ru*r1;
vv = rv*r1;
wv = rw*r1;
E = rE*r1;
q = 0.5*(uv*uv+vv*vv+wv*wv);
p = gam1*(rE-r*q);
p = pmin + lmax(p-pmin,alpha);
dp = atan(alpha*(p - pmin))/pi + (alpha*(p - pmin))/(pi*(alpha^2*(p - pmin)^2 + 1)) + 1/2;

h = E+p*r1;
fi = [ru, ru*uv+p, rv*uv, rw*uv, ru*h, ...
      rv, ru*vv, rv*vv+p, rw*vv, rv*h, ...
      rw, ru*wv, rv*wv, rw*wv+p, rw*h];

ux = (rux - rx*uv)*r1;
vx = (rvx - rx*vv)*r1;
wx = (rwx - rx*wv)*r1;
qx = uv*ux + vv*vx + wv*wx;
px = gam1*(rEx - rx*q - r*qx);
px = px*dp;
%Tx = gam*M2*(px*r - p*rx)*r1^2;
Tx = 1/gam1*(px*r - p*rx)*r1^2;

uy = (ruy - ry*uv)*r1;
vy = (rvy - ry*vv)*r1;
wy = (rwy - ry*wv)*r1;
qy = uv*uy + vv*vy + wv*wy;
py = gam1*(rEy - ry*q - r*qy);
py = py*dp;
%Ty = gam*M2*(py*r - p*ry)*r1^2;
Ty = 1/gam1*(py*r - p*ry)*r1^2;

uz = (ruz - rz*uv)*r1;
vz = (rvz - rz*vv)*r1;
wz = (rwz - rz*wv)*r1;
qz = uv*uz + vv*vz + wv*wz;
pz = gam1*(rEz - rz*q - r*qz);
pz = pz*dp;
%Tz = gam*M2*(pz*r - p*rz)*r1^2;
Tz = 1/gam1*(pz*r - p*rz)*r1^2;

T = p/(gam1*r);
T = Tref/Tinf * T;
mu = getViscosity(muRef,Tref,T,1);
mu = mu + avs;
fc = mu*gam/(Pr);

txx = mu*c23*(2*ux - vy - wz) + (avb+avr)*(ux+vy+wz);
txy = mu*(uy + vx);
txz = mu*(uz + wx);
tyy = mu*c23*(2*vy - ux - wz) + (avb+avr)*(ux+vy+wz);
tyz = mu*(vz + wy);
tzz = mu*c23*(2*wz - ux - vy) + (avb+avr)*(ux+vy+wz);

fv = [0, txx, txy, txz, uv*txx + vv*txy + wv*txz + (fc+avk)*Tx, ...
      0, txy, tyy, tyz, uv*txy + vv*tyy + wv*tyz + (fc+avk)*Ty,...
      0, txz, tyz, tzz, uv*txz + vv*tyz + wv*tzz + (fc+avk)*Tz];
  
f = fi+fv;         

