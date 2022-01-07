function f = getfluxav_3d(udg,odg,param)
gam = param(1);
gam1 = gam - 1.0;
Re = param(2);
Pr = param(3);
Minf = param(4);
Re1 = 1/Re;
M2 = Minf^2;
c23 = 2.0/3.0;
fc = 1/(gam1*M2*Re*Pr);
pmin = 0.2;
rmin = 0.2;
alpha = 1e3;

avb  = odg(1);
avr  = odg(2);
avs  = odg(3);
av   = avb+avr+avs;
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

Tx = gam*M2*(px*r - p*rx)*r1^2;
uy = (ruy - ry*uv)*r1;
vy = (rvy - ry*vv)*r1;
wy = (rwy - ry*wv)*r1;
qy = uv*uy + vv*vy + wv*wy;
py = gam1*(rEy - ry*q - r*qy);
py = py*dp;

Ty = gam*M2*(py*r - p*ry)*r1^2;
uz = (ruz - rz*uv)*r1;
vz = (rvz - rz*vv)*r1;
wz = (rwz - rz*wv)*r1;
qz = uv*uz + vv*vz + wv*wz;
pz = gam1*(rEz - rz*q - r*qz);
pz = pz*dp;

Tz = gam*M2*(pz*r - p*rz)*r1^2;
txx = Re1*c23*(2*ux - vy - wz);
txy = Re1*(uy + vx);
txz = Re1*(uz + wx);
tyy = Re1*c23*(2*vy - ux - wz);
tyz = Re1*(vz + wy);
tzz = Re1*c23*(2*wz - ux - vy);
fv = [0, txx, txy, txz, uv*txx + vv*txy + wv*txz + fc*Tx, ...
      0, txy, tyy, tyz, uv*txy + vv*tyy + wv*tyz + fc*Ty,...
      0, txz, tyz, tzz, uv*txz + vv*tyz + wv*tzz + fc*Tz];
  
fav = [av*rx, av*rux, av*rvx, av*rwx, av*rEx, ...
       av*ry, av*ruy, av*rvy, av*rwy, av*rEy, ...
       av*rz, av*ruz, av*rvz, av*rwz, av*rEz];

f = fi+fv+fav;         
  
