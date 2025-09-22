function [p,t]=uniref(p,t,nref)
%UNIREF 2-D Uniform Mesh Refiner. Subdivides a mesh by half NREF times
%   [P,T]=UNIREF(P,T,NREF)
%
%      P:         Node positions (NP,2)
%      T:         Triangle indices (NT,3)
%      NREF:      Number of Refinements
%

if nargin<3, nref=1; end

    
for iref=1:nref
    np=size(p,1);
    nt=size(t,1);
    pair=[t(:,[1,2]);t(:,[1,3]);t(:,[2,3])];
    [pair,pairi,pairj]=unique(sort(pair,2),'rows');
    pmid=(p(pair(:,1),:)+p(pair(:,2),:))/2;
    t1=t(:,1);
    t2=t(:,2);
    t3=t(:,3);
    t12=pairj(1:nt)+np;
    t13=pairj(nt+1:2*nt)+np;
    t23=pairj(2*nt+1:3*nt)+np;

    t=[t1,t12,t13;
     t12,t23,t13;
     t2,t23,t12;
     t3,t13,t23];
    p=[p;pmid];        
end
