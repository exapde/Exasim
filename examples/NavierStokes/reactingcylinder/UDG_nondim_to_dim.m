function [UDG_dim, UH_dim] = UDG_nondim_to_dim(param, UDG, UH)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
rho_scale   = param{1};
u_scale     = param{2};
rhoe_scale  = param{3};
T_scale     = param{4};
mu_scale    = param{5};
kappa_scale = param{6};
cp_scale    = param{7};
L_scale     = param{8};
Ec = param{9};
% Ec = 1.0;
Pr = param{10};
Re = param{11};
ns  = 5;
UDG_dim = 0*UDG;
UH_dim = 0*UH;
ncu = 8;

UDG_dim(:,1:ns,:) = UDG(:,1:ns,:)*rho_scale;
UDG_dim(:,ncu + (1:ns),:) = UDG(:,ncu + (1:ns),:)*rho_scale / L_scale;
UDG_dim(:,2*ncu + (1:ns),:) = UDG(:,2*ncu + (1:ns),:)*rho_scale / L_scale;

UH_dim(1:ns,:) = UH(1:ns,:)*rho_scale;

UDG_dim(:,ns+1:ns+2,:) = UDG(:,ns+1:ns+2,:) *rho_scale*u_scale;
UDG_dim(:,ncu + (ns+1:ns+2),:) = UDG(:,ncu + (ns+1:ns+2),:)*rho_scale*u_scale / L_scale;
UDG_dim(:,2*ncu + (ns+1:ns+2),:) = UDG(:,2*ncu + (ns+1:ns+2),:)*rho_scale*u_scale / L_scale;

UH_dim(ns+1:ns+2,:) = UH(ns+1:ns+2,:)*rho_scale*u_scale;

UDG_dim(:,ns+3,:) = UDG(:,ns+3,:)*rhoe_scale;
UDG_dim(:,ncu + (ns+3),:) = UDG(:,ncu + (ns+3),:)*rhoe_scale / L_scale;
UDG_dim(:,2*ncu + (ns+3),:) = UDG(:,2*ncu + (ns+3),:)*rhoe_scale / L_scale;

UH_dim(ns+3,:) = UH(ns+3,:)* rhoe_scale;



end