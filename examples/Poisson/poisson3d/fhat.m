function [fh,fh_udg,fh_uh] = fhat(nl,p,udg,uh,param,time,signe)
%FHAT flux function
%   [fh,fhu,fhq,fhm] = fhat(nl,p,u,q,m,param)
%
%      NL(N,ND)              Normal N points
%      P(N,ND)               Coordinates for N points
%      U(N,NC)               Unknown vector for N points with NC components
%      Q(N,NC,ND)            Flux vector for N points with NC components in the
%                            coordinate directions
%      M(N,NC)               Hybrid unkowns for N points with NC components
%      PARAM                 Parameter list
%      FH(N,NC):              Volume flux at N points
%      FHU(N,NC,NC):         Jacobian of the flux flux vector w.r.t. U
%      FHQ(N,NC,NC,ND):      Jacobian of the flux flux vector w.r.t. Q
%      FHQ(N,NC,NC):         Jacobian of the flux flux vector w.r.t. Q

[ng,nc] = size(udg);
nch = 1;
nq = nc-nch;
nd = nq;

kappa = 1;
u     = udg(:,nch);
q     = reshape(udg(:,nch+1:nc),[ng nch nd]);

tau = ones(ng,1);
fh = kappa*reshape(mapContractK(nl,q,[],2,1,3,2,1),[ng nch]) + tau.*(u-uh);

fh_uh = -tau;
fh_u = tau;            
fh_q = kappa*reshape(nl,[ng,nch,nd]);
fh_udg = cat(3,fh_u,fh_q);


