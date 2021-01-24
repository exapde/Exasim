function [p1,t1,c1,p2,t2,c2,p3,c3,t3] = isosurfaceplot(mesh, udg, isoval, cdg, nref)

if mesh.nd ~= 3
    error('isosurface plot is only for 3D field.');
end
if (mesh.porder == 0) || (size(udg,1) == 1)
    error('porder must be greater than 0.');
end

if nargin <= 4
    nref = [];   
end

npl = size(udg,1);
udg = reshape(udg,npl,[]);
cdg = reshape(cdg,npl,[]);

umin = min(udg,[],1);
umax = max(udg,[],1);

x = reshape(mesh.dgnodes(:,1,:),npl,[]);
y = reshape(mesh.dgnodes(:,2,:),npl,[]);
x = mean(x,1);
y = mean(y,1);

% elements contain the isosurface
inde = find((umin<isoval) & (isoval < umax) & (y > 0) & (y < 0.1) & (x > 0.37) & (x<0.65));
length(inde)
pause

if isempty(inde)
    warning('No isosurface exists');
    return;
end

nel = length(inde);

% get scalar fields and coordinates on those elements
udg = udg(:,inde);
cdg = cdg(:,inde);
xdg = reshape(mesh.dgnodes(:,1,inde),[npl nel]);
ydg = reshape(mesh.dgnodes(:,2,inde),[npl nel]);
zdg = reshape(mesh.dgnodes(:,3,inde),[npl nel]);

xmin = min(xdg,[],1);
inde = find(xmin>=0);
udg = udg(:,inde);
cdg = cdg(:,inde);
xdg = xdg(:,inde);
ydg = ydg(:,inde);
zdg = zdg(:,inde);

if mesh.porder>1
    % subdivision
    size(udg)
    [xdg,ydg,zdg,udg,cdg] = scalarrefine(mesh,xdg,ydg,zdg,udg,cdg,nref);
    size(udg)
    
    umin = min(udg,[],1);
    umax = max(udg,[],1);

    % elements contain the isosurface
    inde = find((umin<isoval) & (isoval < umax));
    nel = length(inde);
    npl = size(udg,1);

    % get scalar fields and coordinates on those elements
    udg = udg(:,inde);
    cdg = cdg(:,inde);
    xdg = xdg(:,inde);
    ydg = ydg(:,inde);
    zdg = zdg(:,inde);
end

% if hexes then convert them into tets
if npl==8
    % cube to tets
    c2t = [1 2 4 6; 1 5 4 6; 5 8 4 6; 2 3 4 6; 7 3 4 6; 7 8 4 6]';

    udg = reshape(udg(c2t,:),[4 6*nel]);
    cdg = reshape(cdg(c2t,:),[4 6*nel]);
    xdg = reshape(xdg(c2t,:),[4 6*nel]);
    ydg = reshape(ydg(c2t,:),[4 6*nel]);
    zdg = reshape(zdg(c2t,:),[4 6*nel]);
end

% compute and plot the isosurface
[p1,t1,c1,p2,t2,c2,p3,c3,t3] = marchingtet(xdg, ydg, zdg, udg, isoval, cdg);

function [xref,yref,zref,uref,cref] = scalarrefine(mesh,x,y,z,u,c,nref)

nd = 3;
nt=size(x,2);
porder=mesh.porder;
plocal=mesh.plocal;
tlocal=mesh.tlocal;

if isempty(nref), nref=ceil(log2(max(porder,1))); end
if size(tlocal,2)==4  
    A0=koornwinder(plocal(:,1:nd),porder);
    [plocal,tlocal]=uniref3d(plocal,tlocal,nref);    
    A=koornwinder(plocal(:,1:nd),porder)/A0;
else
    A0=tensorproduct(plocal(:,1:nd),porder);
    %m = porder*(nref+1)+1; 
    m = porder+1+nref; 
    [plocal,tlocal]=cubemesh(m,m,m,1);
    A=tensorproduct(plocal(:,1:nd),porder)/A0;  
end

npln=size(plocal,1);
t = kron(ones(nt,1),tlocal)+kron(npln*(0:nt-1)',0*tlocal+1);
ne = size(t,1);
np = npln*nt;

uref = A*u;
uref = reshape(uref,[np 1]);
uref = reshape(uref(t'), [size(t,2) ne]);

cref = A*c;
cref = reshape(cref,[np 1]);
cref = reshape(cref(t'), [size(t,2) ne]);

xref = A*x;
xref = reshape(xref,[np 1]);
xref = reshape(xref(t'), [size(t,2) ne]);

yref = A*y;
yref = reshape(yref,[np 1]);
yref = reshape(yref(t'), [size(t,2) ne]);

zref = A*z;
zref = reshape(zref,[np 1]);
zref = reshape(zref(t'), [size(t,2) ne]);

% sz=size(mesh.dgnodes(:,1:nd,:)); if length(sz)==2, sz = [sz,1]; end
% dgref=reshape(A*reshape(mesh.dgnodes(:,1:nd,:),npl,sz(2)*sz(3)),[npln,sz(2),sz(3)]);
% dgref = reshape(permute(dgref,[1,3,2]),[np,nd]);
% dgref = reshape(dgref(t',:),[size(t,2) ne nd]);
% 




