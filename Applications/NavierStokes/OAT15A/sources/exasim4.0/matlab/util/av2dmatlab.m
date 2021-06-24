function [avb, avk, avm] = av2dmatlab(udg,odg,param,porder)

gam  = param(1);
gam1 = gam - 1.0;        
Minf = param(4);

% slope regulaization parameter 
alpha = 1.0e2;

% regulaization parameters for the bulk viscosity
bbeta = 0.01;
bmin = 0.0;
bmax = 2.0 * (2.0 / sqrt(gam*gam - 1.0));  
bkh = 1.5;

% regulaization parameters for the thermal viscosity
kbeta = 0.5;
kmin = 0.0;
kmax = 2.0;
kkh = 1.0;

% regulaization parameters for the shear viscosity
mbeta = 0.5;
mmin = 0.0;
mmax = 2.0;
mkh = 1.0;

%
rmin = 1e-6;
pmin = 1e-6;
rEmin = 1e-6;

M0 = odg(:,1,:);
M1 = odg(:,2,:);
M2 = odg(:,3,:);
M3 = odg(:,4,:);
Minv0 = odg(:,5,:);
Minv1 = odg(:,6,:);
Minv2 = odg(:,7,:);
Minv3 = odg(:,8,:);
mh   = odg(:,9,:);

r    = udg(:,1,:);
ru   = udg(:,2,:);
rv   = udg(:,3,:);
rE   = udg(:,4,:);

rx   = udg(:,5,:);
rux  = udg(:,6,:);
rvx  = udg(:,7,:);
rEx  = udg(:,8,:);

ry   = udg(:,9,:);
ruy  = udg(:,10,:);
rvy  = udg(:,11,:);
rEy  = udg(:,12,:);

r = rmin + lmax(r-rmin,alpha);
% dr = atan(alpha.*(r - rmin))./pi + (alpha.*(r - rmin))./(pi.*(alpha.^2.*(r - rmin).^2 + 1)) + 1./2;
% rx = rx.*dr;
% ry = ry.*dr;

rE = rEmin + lmax(rE-rEmin,alpha);
% drE = atan(alpha.*(rE - rEmin))./pi + (alpha.*(rE - rEmin))./(pi.*(alpha.^2.*(rE - rEmin).^2 + 1)) + 1./2;
% rEx = rEx.*drE;
% rEy = rEy.*drE;

r1   = 1./r;
uv   = ru.*r1;
vv   = rv.*r1;
E    = rE.*r1;
q    = 0.5.*(uv.*uv+vv.*vv);
p    = gam1.*(rE-r.*q);
p    = pmin + lmax(p-pmin,alpha);
H    = E+p.*r1;
T    = p ./(gam1.*r);
T0   = T + q./gam;
c_star = sqrt( (2.*gam1.*H) ./ (gam+1) ); 
c_star_infty = (1 ./ Minf) .* sqrt( (2./(gam+1)) .* (1+0.5.*gam1.*Minf.*Minf) );
v_ref = sqrt(1 + c_star_infty.*c_star_infty);

%dp = atan(alpha.*(p - pmin))./pi + (alpha.*(p - pmin))./(pi.*(alpha.^2.*(p - pmin).^2 + 1)) + 1./2;

ux  = (rux - rx.*uv).*r1;
vx  = (rvx - rx.*vv).*r1;    
qx  = uv.*ux + vv.*vx;
px  = gam1.*(rEx - rx.*q - r.*qx);
%px  = px.*dp;
Tx  = (px.*r - p.*rx)./(gam1.*r.*r);

uy  = (ruy - ry.*uv).*r1;
vy  = (rvy - ry.*vv).*r1;    
qy  = uv.*uy + vv.*vy;
py  = gam1.*(rEy - ry.*q - r.*qy);
%py  = py.*dp;
Ty  = (py.*r - p.*ry)./(gam1.*r.*r);
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

nx = Tx;
ny = Ty;
nNorm = sqrt(nx.*nx + ny.*ny + 1.0e-20);
nx = nx ./ nNorm;
ny = ny ./ nNorm;
kh = 1.0./sqrt(Minv0.*nx.*nx + Minv1.*nx.*ny + Minv2.*ny.*nx + Minv3.*ny.*ny + 1.0e-20);
norm_gradT_master = sqrt(M0.*Tx.*Tx + M1.*Tx.*Ty + M2.*Ty.*Tx + M3.*Ty.*Ty);

% thermal viscosity
x = (kkh.*kh./ (porder)) .* norm_gradT_master ./ T0;
f = limiting(x,kmin,kmax,alpha,kbeta);
avk =  r.*gam.*(kkh.*kh./(porder)) .* sqrt(uv.*uv + vv.*vv + c_star.*c_star) .* f;

mh=1;
% shear viscosity
shearOpInNumer = (mh./ (porder)) .* sqrt(S_xx.*S_xx + S_xy.*S_xy + S_yx.*S_yx + S_yy.*S_yy);
x = mkh .* shearOpInNumer ./ v_ref;
f = limiting(x,mmin,mmax,alpha,mbeta);
avm = r.*(mkh.*mh./(porder)) .* sqrt(uv.*uv + vv.*vv + c_star.*c_star) .* f;

% ShearRatio = div_v.*div_v ./ (div_v.*div_v + S_xx.*S_xx + S_xy.*S_xy + S_yx.*S_yx + S_yy.*S_yy + 1.0e-16);
% x = -div_v.*DucrosRatio.*ShearRatio;
% f = limiting(x,bmin,bmax,alpha,bbeta);
% avb = f;

