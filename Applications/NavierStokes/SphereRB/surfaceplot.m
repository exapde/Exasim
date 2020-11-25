function surfaceplot(mesh,master,nref,UDG,clim)

p = mesh.p;
t = mesh.t;
f = mesh.f;

order = master.porder;
npe = master.npe;
perm = master.perm;
ne = size(t,2);
[nd, np] = size(p);
[npf,nfe] = size(perm);

dgnodesm = createdgnodes(p,t,mesh.f,[],[],order);

elcon = zeros(npf,nfe,ne);
for ie=1:ne
    elcon(:,:,ie) = perm + (ie-1)*npe;
end



n = [0 0 1];
p0 =[0, 0, 0];
D  = n*p0';
vis = (abs(n*p - D) < 1.e-4);

face = getelemface(3,1);
f2t = mkf2e(t,1,3);
nf = size(f2t,2);

f = zeros(4,nf);
for i = 1:nf
    f(:,i) = t(face(:,f2t(2,i)),f2t(1,i));
end
    
ib = prod(vis(f),1)>0;
ind = 1:nf;
bf = ind(ib);
ns = length(bf);

dgnodes = zeros(npf,nd,ns);

facecon = faceconnectivity2(t,f2t,3,1,order);


sz = size(UDG);
if length(sz)==2
    u = reshape(UDG(1,:,:),[npf,nf]);
    u = u(:,in);
else
    u = zeros(npf,ns);
end


% shapfc = mkshape(porder,mesh.plocfc,mesh.plocfc,mesh.elemtype);
% dshapft = reshape(permute(shapfc(:,:,2:end),[2 3 1]),[npf*(nd-1) npf]);
j = 1;
for i = 1:nf
    fi = mesh.f(i,end-1:end); % obtain two elements sharing the same face i  
%     if fi(2)>0
%         kf = mesh.t2f(fi(1),:);    % obtain neighboring faces 
%         i1 = find(kf(1,:)==i);  % obtain the index of face i in the 1st element        
%         j1 = elcon(:,i1,fi(1)) - (i-1)*npf;                
%         kf = mesh.t2f(fi(2),:);    % obtain neighboring faces 
%         i2 = find(kf(1,:)==i);  % obtain the index of face i in the second element        
%         j2 = elcon(:,i2,fi(2)) - (i-1)*npf;                
%         [mesh.dgnodes(perm(j1,i1),1:nd,fi(1))-mesh.dgnodes(perm(j2,i2),1:nd,fi(2))]  mesh.tprd                     
%         pause
%     end
    if ismember(fi(2),ib) % face i is a boundary face
        kf = mesh.t2f(fi(1),:);    % obtain neighboring faces 
        i1 = find(kf(1,:)==i);  % obtain the index of face i in the 1st element        
        j1 = elcon(:,i1,fi(1)) - (i-1)*npf;            
        p = mesh.dgnodes(perm(j1,i1),1:nd,fi(1));        
        dgnodes(:,:,j) = p;
        u(:,j) = UDG(perm(j1,i1),1,fi(1));                
        j = j + 1;                
%         dpg = permute(reshape(dshapft*p,[npf nd-1 nd]),[1 3 2]);     
%         nx = reshape(dpg(:,2,1,:).*dpg(:,3,2,:) - dpg(:,3,1,:).*dpg(:,2,2,:),[npf 1]);
%         ny = reshape(dpg(:,3,1,:).*dpg(:,1,2,:) - dpg(:,1,1,:).*dpg(:,3,2,:),[npf 1]);
%         nz = reshape(dpg(:,1,1,:).*dpg(:,2,2,:) - dpg(:,2,1,:).*dpg(:,1,2,:),[npf 1]);            
%         detf = sqrt(nx.^2+ny.^2+nz.^2);
%         nx = nx./detf;
%         ny = ny./detf;
%         nz = nz./detf;                        
%         N = [nx ny nz];
%         U = UDG(perm(j1,i1),:,fi(1));             
    end
end

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

nt = size(u,2);
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


function e1=segcollect(e)
%SEGCOLLECT Collect polygons from edge segments.

ue=unique(e(:));
he=histc(e(:),ue);
current=ue(min(find(he==1))); % Find an endpoint
if isempty(current) % Closed curve
  current=e(1,1);
end
e1=current;
while ~isempty(e)
  ix=min(find(e(:,1)==e1(end)));
  if isempty(ix)
    ix=min(find(e(:,2)==e1(end)));
    if isempty(ix) % >1 disjoint curves, recur
      rest=segcollect(e);
      e1={e1,rest{:}};
      return;
    end
    next=e(ix,1);
  else
    next=e(ix,2);
  end
  e1=[e1,next];
  e(ix,:)=[];
end
e1={e1};




