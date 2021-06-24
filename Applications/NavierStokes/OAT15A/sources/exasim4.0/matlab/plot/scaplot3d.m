function scaplot3d(mesh,u,nref)
%SCAPLOT  Plot Scalar function
%    SCAPLOT(MESH,U,CLIM,NREF,PLTMESH,SURF)
%
%    MESH:       Mesh structure
%    U(NPL,NT):  Scalar fucntion to be plotted
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
  [plocal,tlocal]=uniformlocalpnts(1);
  dgnodes=zeros(4,3,nt);
  dgnodes(:,1,:)=reshape(mesh.p(mesh.t',1),3,1,nt);
  dgnodes(:,2,:)=reshape(mesh.p(mesh.t',2),3,1,nt);
  dgnodes(:,3,:)=reshape(mesh.p(mesh.t',3),3,1,nt);
  u=repmat(u,[4,1]);
else
  plocal=mesh.plocal;
  tlocal=mesh.tlocal;
  dgnodes=mesh.dgnodes;
end
dgnodes = dgnodes(:,1:nd,:);

if isempty(nref), nref=ceil(log2(max(porder,1))); end

if nref>0
  if size(tlocal,2)==4  
    A0=koornwinder(plocal(:,1:nd),porder);
    [plocal,tlocal]=uniref3d(plocal,tlocal,nref);    
    A=koornwinder(plocal(:,1:nd),porder)/A0;
  else
    A0=tensorproduct(plocal(:,1:nd),porder);
    m = porder*(nref+1)+1;
    [plocal,tlocal]=cubemesh(m,m,m,1);
    A=tensorproduct(plocal(:,1:nd),porder)/A0;  
  end
  npln=size(plocal,1);
  sz=size(dgnodes); if length(sz)==2, sz = [sz,1]; end
  dgnodes=reshape(A*reshape(dgnodes,npl,sz(2)*sz(3)),[npln,sz(2),sz(3)]);
  u=A*u;
end

npln=size(plocal,1);
dgnodes=reshape(permute(dgnodes,[1,3,2]),[npln*nt,nd]);
tvis=kron(ones(nt,1),tlocal)+kron(npln*(0:nt-1)',0*tlocal+1);

ne = size(tvis,1);
np = npln*nt;
dgnodes  = [1:np; dgnodes(:,1)'; dgnodes(:,2)'; dgnodes(:,3)'];
tvis  = [(1:ne); tvis'];

fname = sprintf('tmp.inp');
fileID = fopen(fname,'w');
% Header
fprintf(fileID,'%d %d 1 0 0\n',[np,ne]);
% Point coordinates
fprintf(fileID,'%d %1.3f %1.3f %1.3f\n',dgnodes);
% Element connectivities
if mesh.elemtype==0
    fprintf(fileID,'%d 0 tet %d %d %d %d\n',tvis);
else
    fprintf(fileID,'%d 0 hex %d %d %d %d %d %d %d %d\n',tvis);    
end
% Point field
fprintf(fileID,'%d 1\n',1);    
fprintf(fileID,'U,none\n');
%fprintf(fileID,'Uimag,none\n');
fprintf(fileID,'%d %1.3f\n',[1:np; u(:)']);
fclose(fileID);

