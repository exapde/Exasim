function plotudg(dgnodes,u,plocal,tlocal,porder,clim,nref)
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

if nargin<6
    clim = [];
end
if nargin<7
    nref=0;
end

[npl,~,nt]=size(dgnodes);
u=reshape(u,npl,nt);

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
nodesvis=reshape(permute(dgnodes,[1,3,2]),npln*nt,[]);
tvis=kron(ones(nt,1),tlocal)+kron(npln*(0:nt-1)',0*tlocal+1);

patch('vertices',nodesvis,'faces',tvis,'cdata',u, ...
           'facecol','interp','edgec','none');
      
if ~isempty(clim)
  set(gca,'clim',clim);
end
% set(gcf,'rend','z');
% colorbar;
% axis equal;
% drawnow;
