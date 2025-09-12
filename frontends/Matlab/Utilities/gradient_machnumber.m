function [M, Mgrad, Mx, My] = gradient_machnumber(u, q, gam)
   % pde.physicsmu = [gam Re Pr Minf rinf ruinf rvinf rEinf Tref avk avs];

gam1 = gam - 1.0;

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
r1 = 1./r;
uv = ru.*r1;
vv = rv.*r1;
q = 0.5*(uv.*uv+vv.*vv);
p = gam1.*(rE-r.*q);

u2 = sqrt(2*q);
a = sqrt(gam*p./r);
M = u2./a;

ux = (rux - rx.*uv).*r1;
vx = (rvx - rx.*vv).*r1;
qx = uv.*ux + vv.*vx;
u2x = qx./u2;
px = gam1.*(rEx - rx.*q - r.*qx);
ax = gam *(px.*r - rx.*p)./(2*a.*r.^2) ;
Mx = (u2x*a - u2*ax)./(a.^2);

uy = (ruy - ry.*uv).*r1;
vy = (rvy - ry.*vv).*r1;
qy = uv.*uy + vv.*vy;
py = gam1.*(rEy - ry.*q - r.*qy);
u2y = qy./u2;
ay = gam *(py.*r - ry.*p)./(2*a.*r.^2) ;
My = (u2y*a - u2*ay)./(a.^2);

Mgrad = sqrt(Mx.^2 + My.^2); 

end

