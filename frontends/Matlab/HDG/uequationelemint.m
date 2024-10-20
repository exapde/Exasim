function [Ru, BD, pg, udgg, f, f_udg, s, s_udg] = uequationelemint(master,pde,dgnodes,UDG,SH)
% VOLINTND compute volume integrals 

% Obtain dgnodes, Jacobian matrix and determinant at Gauss points
[pg, Xx, jac] = volgeom(master.shapvt,dgnodes);

ne   = size(UDG,2);
nd   = master.nd;
npe  = master.npe;
nge  = master.nge;

nc   = pde.nc;
ncu  = pde.ncu;
arg  = pde.arg;
time = pde.time;
tdep = pde.tdep;
fc_u = pde.fc_u;
source = str2func(pde.source);
flux   = str2func(pde.flux);

% Shap functions and derivatives
shapvt = master.shapvt;
shapvg = reshape(master.shapvg,[npe nge*(nd+1)]);
shapvgdotshapvl = reshape(master.shapvgdotshapvl,[npe*npe nge (nd+1)]);

% DG solution at Gauss points
udgg = reshape(UDG,[npe ne*nc]);
udgg = shapvt(:,:,1)*udgg;
udgg = reshape(udgg,[nge*ne nc]);

% Fluxes and source at Gauss points
[f, f_udg] = flux( pg, udgg, arg, time);
[s, s_udg] = source( pg, udgg, arg, time); 
f     = reshape(f,[nge ne ncu nd]);
f_udg = permute(reshape(f_udg,[nge ne ncu nd nc]),[1 2 3 5 4]); 
s     = reshape(s(:,1:ncu),[nge*ne ncu]);
s_udg = reshape(s_udg(:,1:ncu,1:nc),[nge*ne ncu nc]);

% Update source term for time-dependent problems
if tdep    
    Stn = reshape(SH(:,:,1:ncu),[npe ne*ncu]);
    Stg = shapvt(:,:,1)*Stn;
    Stg = reshape(Stg,[nge*ne ncu]);

    % axis symmetry
    axisflag = isfield(pde,'axisymmetry');
    if axisflag
        xr = pde.axisymmetry.*abs(pg(:,1));
        xr = reshape(xr,[nge*ne 1]);
    else
        xr = ones(nge*ne,1);
    end
    
    % time-derivative coefficients
    dtcoef = isfield(pde,'dtcoef');
    if dtcoef
        dtcoef = pde.dtcoef;
    else
        dtcoef = ones(ncu,1);
    end
    
    for i=1:ncu
        s(:,i) = s(:,i) + dtcoef(i)*xr.*(Stg(:,i)-udgg(:,i)*fc_u);
        s_udg(:,i,i) = s_udg(:,i,i) - dtcoef(i)*fc_u*xr;
    end    
end

% compute wrk and wrl to time with shape functions
wrk = zeros(nge*(nd+1),ne*ncu);
wrl = zeros(nge*(nd+1),ne*ncu*nc);
wrk(1:nge,:) =  reshape(bsxfun(@times,s,jac),[nge ne*ncu]);
wrl(1:nge,:) =  reshape(bsxfun(@times,s_udg,reshape(jac,[nge*ne 1 1])),[nge ne*ncu*nc]);
for i=1:nd
    fk = bsxfun(@times,f(:,:,:,1),Xx(:,:,1,i));
    fl = bsxfun(@times,f_udg(:,:,:,:,1),Xx(:,:,1,i));
    for j=2:nd
        fk = fk + bsxfun(@times,f(:,:,:,j),Xx(:,:,j,i));
        fl = fl + bsxfun(@times,f_udg(:,:,:,:,j),Xx(:,:,j,i));
    end
    wrk(i*nge+1:(i+1)*nge,:) = reshape(fk,[nge ne*ncu]);
    wrl(i*nge+1:(i+1)*nge,:) = reshape(fl,[nge ne*ncu*nc]);
end

% Volume residual
% Xx = [dxi/dx deta/dx; dxi/dy deta/dy]
% [Phi Phi_xi Phi_eta] x [S.*jac; Fx.*Xx(:,:,1,1)+Fy.*Xx(:,:,2,1); Fx.*Xx(:,:,1,2)+Fy.*Xx(:,:,2,2)]
Ru = shapvg*wrk; % [npe nge*(nd+1)] x [nge*(nd+1) ne*ncu] 
Ru = reshape(Ru,[npe ne ncu]); 

% volume matrices
BD = -reshape(shapvgdotshapvl,[npe*npe nge*(nd+1)])*wrl;
BD = reshape(BD,[npe npe ne ncu nc]);

