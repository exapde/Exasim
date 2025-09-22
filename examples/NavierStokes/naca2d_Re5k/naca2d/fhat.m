function [fh,fh_udg,fh_uh] = fhat(nl,p,udg,uh,param,time)
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
%      FHM(N,NC,NC):         Jacobian of the flux flux vector w.r.t. Lambda

[ng,nch] = size(uh);

% if nch==5
%     [fh,fh_udg,fh_uh] = fhat3d(nl,p,udg,uh,param,time);
%     return;
% end

nd = nch-2;    
    
u = udg(:,1:nch);
q = reshape(udg(:,nch+1:end),ng,nch,[]);

% [f,f_udg] = flux(p,udg,param,time);
% 
% %[An,Anm] = getan(nl,uh,param,2);
% An  = zeros(ng,nch,nch);
% for i=1:nch
%     An(:,i,i)=param{end};
% end
% Anm = zeros(ng,nch,nch,nch);
% 
% fh = permute(mapContractK(f,nl,2,3,1,2,[],1) + mapContractK(An,u-uh,2,3,1,2,[],1),[2 1]);
% %fh = permute(mapContractK(An,u-uh,2,3,1,2,[],1),[2 1]);
% 
% fn_udg = mapContractK(f_udg,nl,[2 4],3,1,2,[],1);
% 
% fh_u = permute(fn_udg(:,1:nch,:),[3 1 2])+An;
% fh_q = permute(fn_udg(:,nch+1:3*nch,:),[3 1 2]);
% fh_udg = cat(3,fh_u,fh_q);

% fh_uh = permute(mapContractK(Anm,u-uh,[2 4],3,1,2,[],1),[3 1 2])-An;    

[f,f_uhdg] = flux(p,[uh,udg(:,nch+1:end)],param,time);

%[An,Anm] = getan(nl,uh,param,2);
An  = zeros(ng,nch,nch);
for i=1:nch
    An(:,i,i)=param{end};
end
Anm = zeros(ng,nch,nch,nch);

%fh = permute(mapContractK(f,nl,2,3,1,2,[],1) + mapContractK(An,u-uh,2,3,1,2,[],1),[2 1]);
%fh = permute(mapContractK(An,u-uh,2,3,1,2,[],1),[2 1]);
fh = f(:,:,1).*nl(:,1) + f(:,:,2).*nl(:,2) + param{end}*(u-uh);

fn_udgh = mapContractK(f_uhdg,nl,[2 4],3,1,2,[],1);

fh_u = An;
fh_q = permute(fn_udgh(:,nch+1:(nd+1)*nch,:),[3 1 2]);
fh_udg = cat(3,fh_u,fh_q);

fh_uh = permute(fn_udgh(:,1:nch,:)+mapContractK(Anm,u-uh,[2 4],3,1,2,[],1),[3 1 2])-An;    

