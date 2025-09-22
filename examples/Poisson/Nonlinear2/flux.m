function [f,f_udg] = flux(p,udg,param,time)
%FLUX Volume flux function
%   [f,fu,fq] = flux(p,u,q,param)
%
%      P(N,ND)              Coordinates for N points
%      U(N,NC)              Unknown vector for N points with NC components
%      Q(N,NC,ND)           Flux vector for N points with NC components in the
%                           coordinate directions
%      PARAM                Parameter list
%      F(N,NC,ND):          Volume flux at N points
%      FU(N,NC,ND,NC):      Jacobian of the flux flux vector w.r.t. U
%      FQ(N,NC,ND,NC,ND):   Jacobian of the flux flux vector w.r.t. Q

[ng,nc] = size(udg);
nch = 1;
ncq = nc-nch;
nd = ncq;

q = reshape(udg(:,nch+1:nc),[ng nch nd]);

f = reshape((2 + udg(:,1).*udg(:,1)).*q,[ng nch nd]);

f_u = zeros(ng, nch, nd, nch);
f_u(:,1,1,1) = 2*udg(:,1).*q(:,1,1);
f_u(:,1,2,1) = 2*udg(:,1).*q(:,1,2);

f_q = zeros(ng, nch, nd, ncq);
for d = 1:nd
  f_q(:,:,d,d) = 2 + udg(:,1).*udg(:,1);
end

f_udg = cat(4,f_u,f_q);



