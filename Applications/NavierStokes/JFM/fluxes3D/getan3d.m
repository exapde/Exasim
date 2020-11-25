function An = getan3d(nl,m,gam,cas)
gam1 = gam - 1.0;
nx = nl(1);
ny = nl(2);
nz = nl(3);
r = m(1);
ru = m(2);
rv = m(3);
rw = m(4);
rE = m(5);
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
if cas==0
    Lambda = [ tanh(1e3*(un-a)) , 0 , 0 , 0 , 0 ;...
                     0 , tanh(1e3*(un)) , 0 , 0 , 0 ;...
                     0 , 0 , tanh(1e3*(un)) , 0 , 0 ;...
                     0 , 0 , 0 , tanh(1e3*(un)), 0 ;...
                     0 , 0 , 0 , 0 , tanh(1e3*(un+a)) ];
else
    Lambda = [ abs(un-a) , 0 , 0 , 0 , 0 ;...
                     0 , abs(un) , 0 , 0 , 0 ;...
                     0 , 0 , abs(un) , 0 , 0 ;...
                     0 , 0 , 0 , abs(un) , 0;...
                     0 , 0 , 0 , 0 , abs(un+a) ];
end
E = simplify(K * Lambda * Kinv);
An = simplify(Tinv * E * T);
% function An = getan_3d(nl,m,gam,cas)
%
% gam1 = gam - 1.0;
% nx = nl(1);
% ny = nl(2);
% nz = nl(3);
%
% r = m(1);
% ru = m(2);
% rv = m(3);
% rw = m(4);
% rE = m(5);
%
% r1 = 1/r;
% uv = ru*r1;
% vv = rv*r1;
% wv = rw*r1;
% E = rE*r1;
% q = 0.5*(uv*uv+vv*vv+wv*wv);
% p = gam1*(rE-r*q);
% h = E+p*r1;
% a = sqrt(gam*p*r1);
%
% run = ru*nx + rv*ny + rw*nz;
% un = run/r;
%
% t1x = 0.2039; % Weird direction using irrational numbers to make unlikely it is parallel to normal
% t1y = t1x*pi;
% t1z = t1y*pi/exp(1);
% t1x = t1x - (t1x*nx+t1y*ny+t1z*nz)*nx; % Make t1 orthogonal to normal
% t1y = t1y - (t1x*nx+t1y*ny+t1z*nz)*ny;
% t1z = t1z - (t1x*nx+t1y*ny+t1z*nz)*nz;
% norm_t1 = sqrt(t1x*t1x + t1y*t1y + t1z*t1z);
% t1x = t1x / norm_t1; % Normalize t1
% t1y = t1y / norm_t1;
% t1z = t1z / norm_t1;
%
% t2x = ny*t1z - nz*t1y;
% t2y = nz*t1x - nx*t1z;
% t2z = nx*t1y - ny*t1x;
% norm_t2 = sqrt(t2x*t2x + t2y*t2y + t2z*t2z);
% t2x = t2x / norm_t2;
% t2y = t2y / norm_t2;
% t2z = t2z / norm_t2;
%
% rut1 = ru*t1x + rv*t1y + rw*t1z;
% ut1 = rut1/r;
% rut2 = ru*t2x + rv*t2y + rw*t2z;
% ut2 = rut2/r;
%
%
% K = [ 1 , 1 , 0 , 0, 1 ;...
% un-a , un , 0 , 0 , un+a ;...
% ut1 , ut1 , 1 , 0, ut1 ;...
% ut2 , ut2 , 0 , 1, ut2 ;...
% h - un*a , (1/2)*(un^2 + ut1^2 + ut2^2) , ut1 , ut2 , h+un*a ];
% Kinv = (gam1/(2*a^2))*[ h + (a/gam1)*(un-a) , -(un+a/gam1) , -ut1 , -ut2 , 1 ;...
% -2*h + (4/gam1)*a^2 , 2*un , 2*ut1 , 2*ut2 , -2 ;...
% -2*(ut1*a^2)/gam1 , 0 , 2*(a^2)/gam1 , 0 , 0 ;...
% -2*(ut2*a^2)/gam1 , 0 , 0 , 2*(a^2)/gam1 , 0 ;...
% h - a*(un+a)/gam1 , -un+a/gam1 , -ut1 , -ut2 , 1 ];
%
% T = [ 1 , 0 , 0 , 0 , 0;...
% 0 , nx , ny , nz , 0;...
% 0 , t1x , t1y , t1z , 0;...
% 0 , t2x , t2y , t2z , 0;...
% 0 , 0 , 0 , 0 , 1];
% Tinv = [1, 0, 0, 0, 0;...
% 0, (t1y*t2z - t2y*t1z)/(nx*t1y*t2z - nx*t2y*t1z - ny*t1x*t2z + ny*t2x*t1z + nz*t1x*t2y - nz*t2x*t1y), -(ny*t2z - nz*t2y)/(nx*t1y*t2z - nx*t2y*t1z - ny*t1x*t2z + ny*t2x*t1z + nz*t1x*t2y - nz*t2x*t1y), (ny*t1z - nz*t1y)/(nx*t1y*t2z - nx*t2y*t1z - ny*t1x*t2z + ny*t2x*t1z + nz*t1x*t2y - nz*t2x*t1y), 0;...
% 0, -(t1x*t2z - t2x*t1z)/(nx*t1y*t2z - nx*t2y*t1z - ny*t1x*t2z + ny*t2x*t1z + nz*t1x*t2y - nz*t2x*t1y), (nx*t2z - nz*t2x)/(nx*t1y*t2z - nx*t2y*t1z - ny*t1x*t2z + ny*t2x*t1z + nz*t1x*t2y - nz*t2x*t1y), -(nx*t1z - nz*t1x)/(nx*t1y*t2z - nx*t2y*t1z - ny*t1x*t2z + ny*t2x*t1z + nz*t1x*t2y - nz*t2x*t1y), 0;...
% 0, (t1x*t2y - t2x*t1y)/(nx*t1y*t2z - nx*t2y*t1z - ny*t1x*t2z + ny*t2x*t1z + nz*t1x*t2y - nz*t2x*t1y), -(nx*t2y - ny*t2x)/(nx*t1y*t2z - nx*t2y*t1z - ny*t1x*t2z + ny*t2x*t1z + nz*t1x*t2y - nz*t2x*t1y), (nx*t1y - ny*t1x)/(nx*t1y*t2z - nx*t2y*t1z - ny*t1x*t2z + ny*t2x*t1z + nz*t1x*t2y - nz*t2x*t1y), 0;...
% 0, 0, 0, 0, 1];
%
% if cas==0
% Lambda = [ tanh(1e3*(un-a)) , 0 , 0 , 0 , 0 ;...
% 0 , tanh(1e3*(un)) , 0 , 0 , 0 ;...
% 0 , 0 , tanh(1e3*(un)) , 0 , 0 ;...
% 0 , 0 , 0 , tanh(1e3*(un)), 0 ;...
% 0 , 0 , 0 , 0 , tanh(1e3*(un+a)) ];
% else
% Lambda = [ abs(un-a) , 0 , 0 , 0 , 0 ;...
% 0 , abs(un) , 0 , 0 , 0 ;...
% 0 , 0 , abs(un) , 0 , 0 ;...
% 0 , 0 , 0 , abs(un) , 0;...
% 0 , 0 , 0 , 0 , abs(un+a) ];
% end
% E = (K * Lambda * Kinv);
% An = (Tinv * E * T);
%E = simplify(K * Lambda * Kinv);
%An = simplify(Tinv * E * T);
% gam1 = gam - 1.0;
% nc = length(m);
%
% nx = nl(1);
% ny = nl(2);
%
% r = m(1);
% ru = m(2);
% rv = m(3);
% rE = m(4);
%
% r1 = 1/r;
% uv = ru*r1;
% vv = rv*r1;
% E = rE*r1;
% af = 0.5*(uv*uv+vv*vv);
% p = gam1*(rE -r*af);
% h = E + p*r1;
% c2 = gam* p*r1;
% c = sqrt(c2);
% un = uv*nx + vv*ny;
%
% if absolute
% rlam1 = abs(un+c);
% rlam2 = abs(un-c);
% rlam3 = abs(un);
% else
% rlam1 = un+c;
% rlam2 = un-c;
% rlam3 = abs(un);
% end
%
% s1 = 0.5*(rlam1 + rlam2);
% s2 = 0.5*(rlam1 - rlam2);
% %An = syms(nc,nc);
%
% cc1 = gam1*(s1-rlam3)*af/c2-(s2*un/c);
% cc2 = gam1*s2*af/c-(s1-rlam3)*un;
% An(:,1) = [rlam3+cc1; cc1*uv+cc2*nx; cc1*vv+cc2*ny; cc1*h+cc2*un];
%
% cc1 = -gam1*(s1-rlam3)*uv/c2+(s2*nx/c);
% cc2 = -gam1*s2*uv/c + (s1-rlam3)*nx;
% An(:,2) = [cc1; rlam3+cc1*uv+cc2*nx; cc1*vv+cc2*ny; cc1*h+cc2*un];
%
% cc1 = -gam1*(s1-rlam3)*vv/c2+(s2*ny/c);
% cc2 = -gam1*s2*vv/c+(s1-rlam3)*ny;
%
% An(:,3) = [cc1; cc1*uv+cc2*nx; rlam3+cc1*vv+cc2*ny; cc1*h+cc2*un];
%
% cc1 = gam1*(s1-rlam3)/c2;
% cc2 = gam1*s2/c;
% An(:,4) = [cc1; cc1*uv+cc2*nx; cc1*vv+cc2*ny; rlam3+cc1*h+cc2*un];
%
% % Add the viscous stabilization
% % An(1,1) = An(1,1) + 0;
% % An(2,2) = An(2,2) + 1/(Re);
% % An(3,3) = An(3,3) + 1/(Re);
% % An(4,4) = An(4,4) + 1/(gam1*M2*Re*Pr);
%
%
%
