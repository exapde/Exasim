function [hfx, hfy, Tx, Ty, kappa] = heatflux(u, q, mu)
   % pde.physicsmu = [gam Re Pr Minf rinf ruinf rvinf rEinf Tref avk avs];

gam = mu(1);
gam1 = gam - 1.0;
Re = mu(2);
Pr = mu(3);
Minf = mu(4);
Tref = mu(10);
muRef = 1/Re;
Tinf = 1/(gam*gam1*Minf^2);

% regularization mueters
alpha = 1.0e3;
rmin = 1.0e-2;
pmin = 1.0e-3;

r = u(:,1,:);
ru = u(:,2,:);
rv = u(:,3,:);
rE = u(:,4,:);
rx = q(:,1,:);
rux = q(:,2,:);
rvx = q(:,3,:);
rEx = q(:,4,:);
ry = q(:,5,:);
ruy = q(:,6,:);
rvy = q(:,7,:);
rEy = q(:,8,:);

% Regularization of rho (cannot be smaller than rmin)
r = rmin + lmax(r-rmin,alpha);
% Density sensor
dr = atan(alpha*(r - rmin))/pi + (alpha*(r - rmin))./(pi*(alpha^2*(r - rmin).^2 + 1)) + 1/2;
%dr=1;
rx = rx.*dr;
ry = ry.*dr;
r1 = 1/r;
uv = ru.*r1;
vv = rv.*r1;

q = 0.5.*(uv.*uv+vv.*vv);
p = gam1.*(rE-r.*q);
% Regularization of pressure p (cannot be smaller than pmin)
p = pmin + lmax(p-pmin,alpha);
% Pressure sensor
dp = atan(alpha.*(p - pmin))/pi + (alpha.*(p - pmin))./(pi.*(alpha^2.*(p - pmin).^2 + 1)) + 1/2;

ux = (rux - rx.*uv).*r1;
vx = (rvx - rx.*vv).*r1;
qx = uv.*ux + vv.*vx;
px = gam1.*(rEx - rx.*q - r.*qx);
px = px.*dp;
Tx = 1/gam1.*(px.*r - p.*rx).*r1.^2;

uy = (ruy - ry.*uv).*r1;
vy = (rvy - ry.*vv).*r1;
qy = uv.*uy + vv.*vy;
py = gam1.*(rEy - ry.*q - r.*qy);
py = py.*dp;
Ty = 1/gam1.*(py.*r - p.*ry).*r1.^2;

T = p./(gam1.*r);
Tphys = Tref/Tinf .* T;
mu = getViscosity(muRef,Tref,Tphys,1);
kappa = mu.*gam/(Pr);

hfx = (kappa).*Tx;
hfy = (kappa).*Ty;

end
