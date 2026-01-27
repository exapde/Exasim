function [fh] = fhat(nl,p,udg,uh,param,time,wdg,w_flag)
[ng,nch] = size(uh);
nc = size(udg,2);
if w_flag == 1
[f,~] = flux_dim_2(p,[uh,udg(:,nch+1:end)],param,time,wdg);
else
    [f,~] = flux_dim_2(p,[uh,udg(:,nch+1:end)],param,time);
end
% [f,~] = flux_dim_2(p,[udg(:,1:end)],param,time);
tau = param{end};    
% dtau_duh = 0.0;
fh = permute(mapContractK(f,nl,2,3,1,2,[],1),[2 1]) + tau.*(udg(:,1:nch)-uh);
% disp(norm(fh(:)))
% fh_udg = zeros(ng,nch,nc);
% temp = permute(mapContractK(f_uhdg,nl,[2 4],3,1,2,[],1),[3 1 2]);
% fh_udg(:,:,nch+1:end) = temp(:,:,nch+1:end);    
% fh_uh = temp(:,:,1:nch);
% for i=1:nch
%     fh_udg(:,i,i) = tau;
%     fh_uh(:,i,i) = fh_uh(:,i,i) - tau + dtau_duh.*(udg(:,i)-uh(:,i));
% end    
