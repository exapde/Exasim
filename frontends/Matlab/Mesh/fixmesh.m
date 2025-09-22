function [p,t]=fixmesh(p,t)
%FIXMESH  Remove duplicated/unused nodes and fix element orientation.
%   [P,T]=FIXMESH(P,T)

% NOTE: This function only works for 2D simplices (triangles)

% Remove duplicated nodes:
snap=max(max(p,[],1)-min(p,[],1),[],2)*1024*eps;
[foo,ix,jx]=unique(round(p/snap)*snap,'rows');
p=p(ix,:);
t=jx(t);
if size(t,2) == 1, t = t'; end  % This lines ensures the function works for one element

% Remove nodes that are not contained in t:
[pix,ix,jx]=unique(t);
t=reshape(jx,size(t));
p=p(pix,:);

nv = size(t,2); % # vertices
nd = size(p,2); % # dimensions

if ((nd==2) && (nv==3)) || ((nd==3) && (nv==4))
    v = simpvol(p,t);
    flip=v<0;
    t(flip,[1,2])=t(flip,[2,1]);
end


