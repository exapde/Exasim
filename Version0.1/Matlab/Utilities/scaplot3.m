function scaplot3(mesh,u,clim,nref,pltmesh,surf)
%SCAPLOT  Plot Scalar function
%    SCAPLOT(MESH,U,CLIM,NREF,PLTMESH,SURF)
%
%    MESH:       Mesh structure
%    U(NPL,NT):  Scalar function to be plotted
%                NPL = size(mesh.plocal,1)
%                NT = size(mesh.t,1)
%    CLIM:       CLIM(2) Limits for thresholding (default: no limits)
%    NREF:       Number of refinements used for plotting (default=0)
%    PLTMESH:    0 - do not plot mesh
%                1 - plot mesh with straight edges 
%                2 - plot mesh with curved edges (slow)
%    SURF:       0 - Normal 2D view
%                1 - 3D View
%

nd=mesh.nd;
nt=size(mesh.dgnodes,3);
npl=size(mesh.plocal,1);

u=reshape(u,npl,nt);

porder=mesh.porder;

if porder==0
  nref=0;
  %[plocal,tlocal]=uniformlocalpnts(1);
  [plocal,tlocal] = masternodes(1,nd,mesh.elemtype,mesh.nodetype);
  dgnodes=zeros(3,2,nt);
  dgnodes(:,1,:)=reshape(mesh.p(mesh.t',1),3,1,nt);
  dgnodes(:,2,:)=reshape(mesh.p(mesh.t',2),3,1,nt);
  u=repmat(u,[3,1]);
else
  plocal=mesh.plocal;
  tlocal=mesh.tlocal;
  dgnodes=mesh.dgnodes;
end
dgnodes = dgnodes(:,1:nd,:);

if nargin<4 | isempty(nref), nref=ceil(log2(max(porder,1))); end

if nref>0
  if size(tlocal,2)==3  
    A0=koornwinder(plocal(:,1:2),porder);
    [plocal,tlocal]=uniref(plocal,tlocal,nref);
    A=koornwinder(plocal(:,1:2),porder)/A0;
  else
    A0=tensorproduct(plocal(:,1:2),porder);
    m = porder*(nref+1)+1;
    [plocal,tlocal]=squaremesh(m,m,0,1);
    A=tensorproduct(plocal(:,1:2),porder)/A0;  
  end
  npln=size(plocal,1);
  sz=size(dgnodes); if length(sz)==2, sz = [sz,1]; end
  dgnodes=reshape(A*reshape(dgnodes,npl,sz(2)*sz(3)),[npln,sz(2),sz(3)]);
  u=A*u;
end

npln=size(plocal,1);
nodesvis=reshape(permute(dgnodes,[1,3,2]),[npln*nt,2]);
tvis=kron(ones(nt,1),tlocal)+kron(npln*(0:nt-1)',0*tlocal+1);

if nargin>=6 & ~isempty(surf)
   nodesvis = [nodesvis,reshape(u,size(nodesvis(:,1)))];
end

e=boundedges(nodesvis,tvis,mesh.elemtype);
dgnodesltx=nodesvis(:,1);
dgnodeslty=nodesvis(:,2);

nodesvis(:,3) = u(:);

dgnodesltz=u(:);
line(dgnodesltx(e'),dgnodeslty(e'),dgnodesltz(e'),'color',[0,0,0],'LineWidth',1);    
patch('vertices',nodesvis,'faces',tvis,'cdata',u, ...
           'facecol','interp','edgec','none');
if nargin>=3 & ~isempty(clim)
  set(gca,'clim',clim);
end       
view(3);

function e=boundedges(p,t,elemtype)
%BOUNDEDGES Find boundary edges from triangular mesh
%   E=BOUNDEDGES(P,T)

% Form all edges, non-duplicates are boundary edges
if elemtype==0
    edges=[t(:,[1,2]);
           t(:,[1,3]);
           t(:,[2,3])];
    node3=[t(:,3);t(:,2);t(:,1)];
else
    edges=[t(:,[1,2]);
           t(:,[2,3]);
           t(:,[3,4]);
           t(:,[4,1]);];
    node3=[t(:,4);t(:,3);t(:,2);t(:,1)];    
end

edges=sort(edges,2);
[foo,ix,jx]=unique(edges,'rows');
vec=histc(jx,1:max(jx));
qx=find(vec==1);
e=edges(ix(qx),:);
node3=node3(ix(qx));

% Orientation
v1=p(e(:,2),:)-p(e(:,1),:);
v2=p(node3,:)-p(e(:,1),:);
ix=find(v1(:,1).*v2(:,2)-v1(:,2).*v2(:,1)>0);
e(ix,[1,2])=e(ix,[2,1]);
