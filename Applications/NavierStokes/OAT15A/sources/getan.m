function An = getan(nl,m,gam,cas)
gam1 = gam - 1.0;
nx = nl(1);
ny = nl(2);
r = m(1);
ru = m(2);
rv = m(3);
rE = m(4);
r1 = 1/r;
uv = ru*r1;
vv = rv*r1;
E = rE*r1;
q = 0.5*(uv*uv+vv*vv);
p = gam1*(rE-r*q);
h = E+p*r1;
a = sqrt(gam*p*r1);
run = ru*nx + rv*ny;
rut = -ru*ny + rv*nx;
un = run/r;
ut = rut/r;
K = [ 1 , 1 , 0 , 1 ;...
      un-a , un , 0 , un+a ;...
      ut , ut , 1 , ut ;...
      h - un*a , (1/2)*(un^2 + ut^2) , ut , h+un*a ];
Kinv = (gam1/(2*a^2))*[ h + (a/gam1)*(un-a) , -(un+a/gam1) , -ut , 1 ;...
                        -2*h + (4/gam1)*a^2 , 2*un , 2*ut , -2 ;...
                        -2*(ut*a^2)/gam1 , 0 , 2*(a^2)/gam1 , 0 ;...
                        h - a*(un+a)/gam1 , -un+a/gam1 , -ut , 1 ];
T = [ 1 , 0 , 0 , 0;...
      0 , nx , ny , 0;...
      0 , -ny , nx , 0;...
      0 , 0 , 0 , 1];
Tinv = [ 1 , 0 , 0 , 0;...
         0 , nx ,-ny , 0;...
         0 , ny , nx , 0;...
         0 , 0 , 0 , 1];
if cas==0
Lambda = [ tanh(1e3*(un-a)) , 0 , 0 , 0 ;...
                 0 , tanh(1e3*(un)) , 0 , 0 ;...
                 0 , 0 , tanh(1e3*(un)) , 0 ;...
                 0 , 0 , 0 , tanh(1e3*(un+a)) ];
else
Lambda = [ abs(un-a) , 0 , 0 , 0 ;...
                 0 , abs(un) , 0 , 0 ;...
                 0 , 0 , abs(un) , 0 ;...
                 0 , 0 , 0 , abs(un+a) ];
end
E = simplify(K * Lambda * Kinv);
An = simplify(Tinv * E * T);
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
