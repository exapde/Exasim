function [An,An1] = getan3d(n,u)

gam = 1.4;
gam1 = gam - 1.0;
nx = n(1);
ny = n(2);
nz = n(3);
r = u(1);
ru = u(2);
rv = u(3);
rw = u(4);
rE = u(5);
r1 = 1/r;
uv = ru*r1;
vv = rv*r1;
wv = rw*r1;
E = rE*r1;
af = 0.5*(uv*uv+vv*vv+wv*wv);
p = gam1*(rE -r*af);
h = E + p*r1;
c2 = gam*p*r1;
c = sqrt(c2);
un = uv*nx + vv*ny + wv*nz;

rlam1 = tanh(1e3*(un+c));
rlam2 = tanh(1e3*(un-c));
rlam3 = tanh(1e3*(un));

s1 = 0.5*(rlam1 + rlam2);
s2 = 0.5*(rlam1 - rlam2);

An = zeros(5,5);
cc1 = gam1*(s1-rlam3)*af/c2-(s2*un/c);
cc2 = gam1*s2*af/c-(s1-rlam3)*un;
An1 = [rlam3+cc1; cc1*uv+cc2*nx; cc1*vv+cc2*ny; cc1*wv+cc2*nz; cc1*h+cc2*un];
cc1 = -gam1*(s1-rlam3)*uv/c2+(s2*nx/c);
cc2 = -gam1*s2*uv/c + (s1-rlam3)*nx;
An2 = [cc1; rlam3+cc1*uv+cc2*nx; cc1*vv+cc2*ny; cc1*wv+cc2*nz; cc1*h+cc2*un];
cc1 = -gam1*(s1-rlam3)*vv/c2+(s2*ny/c);
cc2 = -gam1*s2*vv/c+(s1-rlam3)*ny;
An3 = [cc1; cc1*uv+cc2*nx; rlam3+cc1*vv+cc2*ny; cc1*wv+cc2*nz; cc1*h+cc2*un];
cc1 = -gam1*(s1-rlam3)*wv/c2+(s2*nz/c);
cc2 = -gam1*s2*wv/c+(s1-rlam3)*nz;
An4 = [cc1; cc1*uv+cc2*nx; cc1*vv+cc2*ny; rlam3+cc1*wv+cc2*nz; cc1*h+cc2*un];
cc1 = gam1*(s1-rlam3)/c2;
cc2 = gam1*s2/c;
An5 = [cc1; cc1*uv+cc2*nx; cc1*vv+cc2*ny; cc1*wv+cc2*nz; rlam3+cc1*h+cc2*un];

An = [An1, An2, An3, An4, An5];

cy = sqrt(nx^2+ny^2);
sy = nz;
cz = nx/cy;
sz = ny/cy;
Tz = [1, 0, 0, 0, 0;...
      0, cz, sz, 0, 0;...
      0, -sz, cz, 0, 0;...
      0, 0, 0, 1, 0;
      0, 0, 0, 0, 1];
Ty = [1, 0, 0, 0, 0;...
      0, cy, 0, sy, 0;...
      0, 0, 1, 0, 0;...
      0, -sy, 0, cy, 0;
      0, 0, 0, 0, 1];
  
T = Ty*Tz;
Tinv = (Tz.')*(Ty.');
U = [r, ru, rv, rw, rE].';
TU = T*U;
run = TU(2);
rvn = TU(3);
rwn = TU(4);
un = run/r;
vn = rvn/r;
wn = rwn/r;
P = gam1*(rE - r*(1/2)*(un^2 + vn^2 + wn^2));
a = sqrt(gam*P/r);
H = rE/r + P/r;
K = [ 1 , 1 , 0 , 0, 1 ;...
      un-a , un , 0 , 0, un+a ;...
      vn , vn , 1 , 0, vn ;...
      wn , wn , 0 , 1, wn ;...
      H - un*a , (1/2)*(un^2+vn^2+wn^2) , vn , wn, H+un*a ];
Kinv = (gam1/(2*a^2))*[ H + (a/gam1)*(un-a) , -(un+a/gam1) , -vn , -wn, 1 ;...
                        -2*H + (4/gam1)*a^2 , 2*un , 2*vn , 2*wn, -2 ;...
                        -2*(vn*a^2)/gam1 , 0 , 2*(a^2)/gam1 , 0, 0 ;...
                        -2*wn*a^2/gam1 , 0 , 0 , 2*a^2/gam1, 0
                        H - a*(un+a)/gam1 , -un+a/gam1 , -vn , -wn, 1 ];

    Lambda = [ tanh(1e3*(un-a)) , 0 , 0 , 0 , 0 ;...
                     0 , tanh(1e3*(un)) , 0 , 0 , 0 ;...
                     0 , 0 , tanh(1e3*(un)) , 0 , 0 ;...
                     0 , 0 , 0 , tanh(1e3*(un)), 0 ;...
                     0 , 0 , 0 , 0 , tanh(1e3*(un+a)) ];


E = K * Lambda * Kinv;
An1 = Tinv * E * T;




