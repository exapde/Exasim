function div = divergence(udg, href)
ns = 5;
nch = 8;
ns = 5;

rho_i = udg(:,1:ns,:);
rho   = sum(rho_i,2);
rhou  = udg(:,ns+1,:);
rhov  = udg(:,ns+2,:);
rhoE  = udg(:,ns+3,:);

drho_dx_i = udg(:,nch + (1:ns),:);
drhou_dx  = udg(:,nch + ns+1,:);
drhov_dx  = udg(:,nch + ns+2,:);
drhoE_dx  = udg(:,nch + ns+2+1,:);
drho_dy_i = udg(:,nch + (nch+1:nch+ns),:);
drhou_dy  = udg(:,nch + nch+ns+1,:);
drhov_dy  = udg(:,nch + nch+ns+2,:);
drhoE_dy  = udg(:,nch + nch+ns+2+1,:);

rhoinv = 1.0 ./ rho;
uv = rhou .* rhoinv; %velocity
vv = rhov .* rhoinv;
E = rhoE .* rhoinv; %energy
uTu2   = 0.5*(uv.*uv+vv.*vv);

% r_i    = udg(:,1:5,:);
% r = sum(r_i,2);

% ru   = udg(:,ns+1,:);
% rv   = udg(:,ns+2,:);
% rx   = udg(:,5,:);
% rux  = udg(:,6,:);
% ry   = udg(:,9,:);    
% rvy  = udg(:,11,:);    

% r1   = 1./r;
% u    = ru.*r1;
% v    = rv.*r1;
ux = (drhou_dx - uv .* sum(drho_dx_i,2)) .* rhoinv;
vy = (drhov_dy - vv .* sum(drho_dy_i,2)) .* rhoinv;
% ux  = (rux - rx.*u).*r1;    
% vy  = (rvy - ry.*v).*r1;            
div = href*(ux+vy);    

