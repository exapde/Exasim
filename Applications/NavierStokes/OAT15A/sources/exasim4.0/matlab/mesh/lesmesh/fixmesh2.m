function [ p,t ] = fixmesh2( p,t )
%FIXMESH2 Remove duplicated/unused nodes (assume element orient correct)
%   [P,T]=FIXMESH2(P,T)

% Remove duplicated nodes:
snap=max(max(p,[],1)-min(p,[],1),[],2)*1024*eps;
[~,ix,jx]=unique(round(p/snap)*snap,'rows');
p=p(ix,:);
t=jx(t);

end

