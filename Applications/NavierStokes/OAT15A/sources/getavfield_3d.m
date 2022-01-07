function avField = getavfield_3d(udg,odg,param,porder)
gam = param(1);
gam1 = gam - 1.0;
Minf = param(4);

rmin = param(5);
pmin = param(6);
avr = param(7);
avs = param(8);
avbulk = param(9);

% slope regulaization parameter
alpha = 1.0e3;

% regulaization parameters for the bulk viscosity
bbeta = 0.25;
bmin = 0.0;
bmax = 2.0 * (2.0 / sqrt(gam*gam - 1.0));
bkh = 1.5;

% regulaization parameters for the shear viscosity
sbeta = 6.0;
smin = 0.0;
smax = 8.0;

avfix = odg(4);
Minv0 = odg(5);
Minv1 = odg(6);
Minv2 = odg(7);
Minv3 = odg(8);
Minv4 = odg(9);
Minv5 = odg(10);
Minv6 = odg(11);
Minv7 = odg(12);
Minv8 = odg(13);
hmin = odg(14);
avfbk = odg(15); 

r = udg(1);
ru = udg(2);
rv = udg(3);
rw = udg(4);
rE = udg(5);
rx = udg(6);
rux = udg(7);
rvx = udg(8);
rwx = udg(9);
ry = udg(11);
ruy = udg(12);
rvy = udg(13);
rwy = udg(14);
rz = udg(16);
ruz = udg(17);
rvz = udg(18);
rwz = udg(19);
r0 = 0.3;
rreg = r0 + lmax(r-r0,alpha);
dr = atan(alpha*(r - r0))/pi + (alpha*(r - r0))/(pi*(alpha^2*(r - r0)^2 + 1)) + 1/2;
%rreg = r;
r1 = 1/rreg;
uv = ru*r1;
vv = rv*r1;
wv = rw*r1;
%E = rE*r1;
q = 0.5*(uv*uv+vv*vv+wv*wv);
%p = gam1*(rE-rreg*q);
%p = pmin + lmax(p-pmin,alpha);
%H = E+p*r1;
%c_star = sqrt( (2*gam1*H) / (gam+1) );
c_star = (1 / Minf) * sqrt( (2/(gam+1)) * (1+0.5*gam1*Minf*Minf) );
ux = (rux - rx*uv)*r1;
vx = (rvx - rx*vv)*r1;
wx = (rwx - rx*wv)*r1;
uy = (ruy - ry*uv)*r1;
vy = (rvy - ry*vv)*r1;
wy = (rwy - ry*wv)*r1;
uz = (ruz - rz*uv)*r1;
vz = (rvz - rz*vv)*r1;
wz = (rwz - rz*wv)*r1;
S_xx = 0.5 * (ux + ux);
S_xy = 0.5 * (uy + vx);
S_xz = 0.5 * (uz + wx);
S_yx = 0.5 * (vx + uy);
S_yy = 0.5 * (vy + vy);
S_yz = 0.5 * (vz + wy);
S_zx = 0.5 * (wx + uz);
S_zy = 0.5 * (wy + vz);
S_zz = 0.5 * (wz + wz);
div_v = - (ux + vy + wz);
vort_x = - (wy - vz);
vort_y = - (uz - wx);
vort_z = - (vx - uy);
vort = sqrt(vort_x*vort_x + vort_y*vort_y + vort_z*vort_z);
shear = sqrt(S_xx*S_xx + S_yx*S_yx + S_zx*S_zx + ...
                     S_xy*S_xy + S_yy*S_yy + S_zy*S_zy + ...
                     S_xz*S_xz + S_yz*S_yz + S_zz*S_zz);

% limit 
sigm = 3.0e3;
div_v = limiting(div_v,-sigm,sigm,alpha,-sigm);
vort = limiting(vort,0.0,sigm,alpha,0.0);
shear = limiting(shear,0.0,sigm,alpha,0.0);

DucrosRatio = div_v*div_v / (div_v*div_v + vort*vort + 1.0e-16);
ShearRatio = div_v*div_v / (div_v*div_v + shear*shear + 1.0e-16);

rx = rx*dr;
ry = ry*dr;
rz = rz*dr;
nx = rx;
ny = ry;
nz = rz;
nNorm = sqrt(nx*nx + ny*ny + nz*nz + 1.0e-16);
nx = nx / nNorm;
ny = ny / nNorm;
nz = nz / nNorm;
bh = 1.0/sqrt(Minv0*nx*nx + Minv1*ny*nx + Minv2*nz*nx + ...
              Minv3*nx*ny + Minv4*ny*ny + Minv5*nz*ny + ...
              Minv6*nx*nz + Minv7*ny*nz + Minv8*nz*nz + 1.0e-16);
bh = limiting(bh,0.0,20.0,alpha,0.0);
          
% bulk viscosity
x = - (bkh*bh/ (porder)) * div_v / c_star;
x = x*DucrosRatio*ShearRatio;
f = limiting(x,bmin,bmax,alpha,bbeta);
avb = (bkh*bh/(porder)) * sqrt(uv*uv + vv*vv + wv*wv + c_star*c_star) * f;
avb = avbulk*avb*avfbk;

% artificial viscosity for continuity equation
p = gam1*(rE-rreg*q);
avd = -avr*(lmin(r-rmin,alpha)+lmin(p-pmin,alpha))*avfbk;
%avd = avd*avfix;

% shear viscosity
x = (0.01/porder)*shear;
f = limiting(x,smin,smax,alpha,sbeta);
f = sqrt(uv*uv + vv*vv + wv*wv + c_star*c_star)*f;
avs = avs*hmin*f*avfix;

avField(1) = avb;
avField(2) = avd;
avField(3) = avs;
