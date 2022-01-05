function solfaceplot(mesh,dgnodes,u,nref,clim)

porder = mesh.porder;
perm = mesh.perm;
nd = mesh.nd;
npf = size(perm,1);
nt = size(dgnodes,3);

plocal=mesh.plocfc;
tlocal=mesh.tlocfc;
if nref>0
  if size(tlocal,2)==3  
    A0=koornwinder(plocal(:,1:2),porder);
    [plocal,tlocal]=uniref(plocal,tlocal,nref);
    A=koornwinder(plocal(:,1:2),porder)/A0;
  else
    A0=tensorproduct(plocal(:,1:2),porder);
    m = porder*(nref+1);
    [plocal,tlocal]=squaremesh(m,m,0,1);
    A=tensorproduct(plocal(:,1:2),porder)/A0;  
  end
  npln=size(plocal,1);
  sz=size(dgnodes); if length(sz)==2, sz = [sz,1]; end
  dgnodes=reshape(A*reshape(dgnodes,npf,sz(2)*sz(3)),[npln,sz(2),sz(3)]);
  u=A*u;
end

npln=size(plocal,1);
nodesvis=reshape(permute(dgnodes,[1,3,2]),[npln*nt,nd]);
tvis=kron(ones(nt,1),tlocal)+kron(npln*(0:nt-1)',0*tlocal+1);

figure(1); clf;
set(axes,'FontSize',16);

patch('vertices',nodesvis,'faces',tvis,'cdata',u(:), ...
           'facecol','interp','edgec','none');
      
if nargin>=5 && ~isempty(clim)
  set(gca,'clim',clim);
end
colorbar('FontSize',15,'Location','East');
axis equal;axis tight;view(3);



