function avb = avbulk2d(udg,odg,param,porder)

gam  = param(1);
gam1 = gam - 1.0;        

% slope regulaization parameter 
alpha = 1.0e3;

% regulaization parameters for the bulk viscosity
bbeta = 0.01;
bmin = 0.0;
bmax = 2.0 * (2.0 / sqrt(gam*gam - 1.0));  
bkh = 1.5;

avfix = odg(:,3,:);
Minv0 = odg(:,4,:);
Minv1 = odg(:,5,:);
Minv2 = odg(:,6,:);
Minv3 = odg(:,7,:);

r    = udg(:,1,:);
ru   = udg(:,2,:);
rv   = udg(:,3,:);
rE   = udg(:,4,:);

rx   = udg(:,5,:);
rux  = udg(:,6,:);
rvx  = udg(:,7,:);

ry   = udg(:,9,:);
ruy  = udg(:,10,:);
rvy  = udg(:,11,:);

r1   = 1./r;
uv   = ru.*r1;
vv   = rv.*r1;
E    = rE.*r1;
q    = 0.5.*(uv.*uv+vv.*vv);
p    = gam1.*(rE-r.*q);
H    = E+p.*r1;
c_star = sqrt( (2.*gam1.*H) ./ (gam+1) ); 
%c_star_infty = (1 ./ Minf) .* sqrt( (2./(gam+1)) .* (1+0.5.*gam1.*Minf.*Minf) );

ux  = (rux - rx.*uv).*r1;
vx  = (rvx - rx.*vv).*r1;    
uy  = (ruy - ry.*uv).*r1;
vy  = (rvy - ry.*vv).*r1;    
S_xx = 0.5 .* (ux + ux);
S_xy = 0.5 .* (uy + vx);
S_yx = 0.5 .* (vx + uy);
S_yy = 0.5 .* (vy + vy);
div_v = - (ux + vy);
vort = - (vx - uy);
DucrosRatio = div_v.*div_v ./ (div_v.*div_v + vort.*vort + 1.0e-16);
ShearRatio = div_v.*div_v ./ (div_v.*div_v + S_xx.*S_xx + S_xy.*S_xy + S_yx.*S_yx + S_yy.*S_yy + 1.0e-16);

nx = rx;
ny = ry;
nNorm = sqrt(nx.*nx + ny.*ny + 1.0e-16);
nx = nx ./ nNorm;
ny = ny ./ nNorm;
bh = 1.0./sqrt(Minv0.*nx.*nx + Minv1.*nx.*ny + Minv2.*ny.*nx + Minv3.*ny.*ny + 1.0e-20);

% bulk viscosity
x = - (bkh.*bh./ (porder)) .* div_v ./ c_star;
x = x.*DucrosRatio.*ShearRatio;
f = limiting(x,bmin,bmax,alpha,bbeta);
avb = r.*(bkh.*bh./(porder)) .* sqrt(uv.*uv + vv.*vv + c_star.*c_star) .* f;
avb = avb.*avfix;

