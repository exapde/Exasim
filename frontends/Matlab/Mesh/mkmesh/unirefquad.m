function [pref,tref]=unirefquad(p,t)
%UNIREF 2-D Uniform Mesh Refiner. Subdivides a mesh by half NREF times
%   [P,T]=UNIREF(P,T,NREF)
%
%      P:         Node positions (NP,2)
%      T:         Triangle indices (NT,3)


nd = 2;
[ne,nv] = size(t);
pe = reshape(p(t',:),[nv ne nd]);
pc = reshape(mean(pe,1),[ne nd]);
pm = 0*pe;
pm(1,:,:) = 0.5*(pe(1,:,:) + pe(2,:,:));
pm(2,:,:) = 0.5*(pe(2,:,:) + pe(3,:,:));
pm(3,:,:) = 0.5*(pe(3,:,:) + pe(4,:,:));
pm(4,:,:) = 0.5*(pe(4,:,:) + pe(1,:,:));

pref = zeros(2*nv+1, ne, nd);
pref(1:nv,:,:) = pe;
pref((nv+1):2*nv,:,:) = pm;
pref(end,:,:) = pc;

% 4-----7-----3
% |     |     |
% |     |     |
% 8-----9-----6
% |     |     | 
% |     |     | 
% 1-----5-----2
tref = zeros(nv, nv, ne);
a = 0:9:9*(ne-1);
tref(1,1,:) = 1 + a;
tref(2,1,:) = 5 + a;
tref(3,1,:) = 9 + a;
tref(4,1,:) = 8 + a;

tref(1,2,:) = 5 + a;
tref(2,2,:) = 2 + a;
tref(3,2,:) = 6 + a;
tref(4,2,:) = 9 + a;

tref(1,3,:) = 9 + a;
tref(2,3,:) = 6 + a;
tref(3,3,:) = 3 + a;
tref(4,3,:) = 7 + a;

tref(1,4,:) = 8 + a;
tref(2,4,:) = 9 + a;
tref(3,4,:) = 7 + a;
tref(4,4,:) = 4 + a;

pref = reshape(pref, [(2*nv+1)*ne, nd]);
tref = reshape(tref, [nv nv*ne]);
tref = tref';

[pref,tref]=fixmesh(pref,tref,1e-6);

figure(1); clf; simpplot(p, t);
figure(2); clf; simpplot(pref, tref);


% % np=size(p,1);
% % nt=size(t,1);
% 
% % p(t',:);
% 
% pair=[t(:,[1,2]);t(:,[1,3]);t(:,[2,3])];
% [pair,pairi,pairj]=unique(sort(pair,2),'rows');
% pmid=(p(pair(:,1),:)+p(pair(:,2),:))/2;
% t1=t(:,1);
% t2=t(:,2);
% t3=t(:,3);
% t12=pairj(1:nt)+np;
% t13=pairj(nt+1:2*nt)+np;
% t23=pairj(2*nt+1:3*nt)+np;
% 
% t=[t1,t12,t13;
%  t12,t23,t13;
%  t2,t23,t12;
%  t3,t13,t23];
% p=[p;pmid];        
