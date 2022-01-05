function faceplot(mesh,nref,UDG,func,clim)

porder = mesh.porder;
perm = mesh.perm;
ne = mesh.ne;
nf = mesh.nf;
nd = mesh.nd;
[npf,nfe] = size(perm);
elcon = reshape(mesh.elcon,[npf nfe ne]);

in = mesh.f(:,end)<0;
ns = sum(in); 
dgnodes = zeros(npf,nd,ns);
u = zeros(npf,ns);

% shapfc = mkshape(porder,mesh.plocfc,mesh.plocfc,mesh.elemtype);
% dshapft = reshape(permute(shapfc(:,:,2:end),[2 3 1]),[npf*(nd-1) npf]);
%UDG = reshape(UDG,[npf nf]);
j = 0;
for i = 1:nf
    fi = mesh.f(i,end-1:end); % obtain two elements sharing the same face i  
    tf = mesh.f(i,1:end-2);
    pf = (mesh.p(tf,:));
    in = feval(func,pf);      
    if in == 1
        j = j + 1;        
        %u(:,j) = UDG(:,i);
        if fi(2)>0
            kf = mesh.t2f(fi(1),:);  % obtain neighboring faces 
            i1 = find(kf(1,:)==i);  % obtain the index of face i in the 1st element        
            j1 = elcon(:,i1,fi(1)) - (i-1)*npf;                
            kf = mesh.t2f(fi(2),:);    % obtain neighboring faces 
            i2 = find(kf(1,:)==i);  % obtain the index of face i in the second element        
            j2 = elcon(:,i2,fi(2)) - (i-1)*npf;                
            u1 = UDG(perm(j1,i1),1,fi(1));
            u2 = UDG(perm(j2,i2),1,fi(2));
            u(:,j) = (u1+u2)/2;        
            dgnodes(:,:,j) = mesh.dgnodes(perm(j1,i1),1:nd,fi(1));                      
            %[mesh.dgnodes(perm(j1,i1),1:nd,fi(1))-mesh.dgnodes(perm(j2,i2),1:nd,fi(2))]                       
            %pause
        else
            kf = mesh.t2f(fi(1),:);    % obtain neighboring faces 
            i1 = find(kf(1,:)==i);  % obtain the index of face i in the 1st element        
            j1 = elcon(:,i1,fi(1)) - (i-1)*npf;            
            u(:,j) = UDG(perm(j1,i1),1,fi(1));
            dgnodes(:,:,j) = mesh.dgnodes(perm(j1,i1),1:nd,fi(1));                                                         
        end
    end
end
dgnodes = dgnodes(:,:,1:j);
dg = dgnodes;
u = u(:,1:j);

plocal=mesh.plocfc;
tlocal=mesh.tlocfc;
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
  dgnodes=reshape(A*reshape(dgnodes,npf,sz(2)*sz(3)),[npln,sz(2),sz(3)]);
  u=A*u;
end



nt = size(u,2);
npln=size(plocal,1);
nodesvis=reshape(permute(dgnodes,[1,3,2]),[npln*nt,nd]);
tvis=kron(ones(nt,1),tlocal)+kron(npln*(0:nt-1)',0*tlocal+1);

%figure(1); clf;
for i=1:3
    x = nodesvis(:,i);
    x = x-mean(x);
    if max(abs(x(:)))<1e-10        
        break;
    end
end
nodesvis(:,i) = [];

% average u
% nn = length(u(:));
% for ii = 1:nn
%     x = nodesvis(ii,:);
%     s = (nodesvis(:,1)-x(1)).^2 + (nodesvis(:,2)-x(2)).^2;
%     ind = find(abs(sqrt(s))<1e-8);
%     if isempty(ind)==0
%         ind
%         u(ii) = mean(u(ind));
%     end
% end

%[nodesvis,tvis]=fixmesh(nodesvis,tvis);

patch('vertices',nodesvis,'faces',tvis,'cdata',u(:), ...
           'facecol','interp','edgec','none');

% patch('vertices',nodesvis,'faces',tvis,'cdata',0*u(:), ...
%            'facecol','interp','edgec','k');
       
if nargin>=5 && ~isempty(clim)
  set(gca,'clim',clim);
end
%view([45 45]);

dgnodes = dg;
dgnodes(:,i,:) = [];
pars={'facecolor',[.8,1,.8],'edgecolor','k','Linew',0.5};
figure(2); clf;
e=boundedges(mesh.plocfc,mesh.tlocfc,mesh.elemtype);
e1=segcollect(e);        
axis equal,axis off
nt=size(dgnodes,3);
hh=zeros(nt,1);
for it=1:nt
    px=dgnodes(:,1,it);
    py=dgnodes(:,2,it);
    pz=0*px;
    hh(it)=patch(px(e1{1}'),py(e1{1}'),pz(e1{1}'),0.0*e1{1}',pars{:});
end        
view(2),axis equal;


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






