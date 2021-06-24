function [pvis,tvis,u] = surfaceplot(mesh,UDG,ib,nref,clim)

porder = mesh.porder;
perm = mesh.perm;
ne = mesh.ne;
nd = mesh.nd;
[npf,nfe] = size(perm);
%elcon = reshape(elcon,[npf nfe ne]);

in = find(mesh.f(:,end)==ib);
if isempty(in)
    error('Boundary is invalid.');
end
ns = length(in); 
dgnodes = zeros(npf,nd,ns);
u = zeros(npf,ns);

for j = 1:ns
    i = in(j);
    fi = mesh.f(i,end-1:end); % obtain two elements sharing the same face i      
    kf = mesh.t2f(fi(1),:);    % obtain neighboring faces 
    i1 = find(kf(1,:)==i);  % obtain the index of face i in the 1st element        
    %j1 = elcon(:,i1,fi(1)) - (i-1)*npf;  
    j1 = 1:1:npf;
    p = mesh.dgnodes(perm(j1,i1),1:nd,fi(1));        
    dgnodes(:,:,j) = p;
    u(:,j) = UDG(perm(j1,i1),1,fi(1));                                 
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
pvis=reshape(permute(dgnodes,[1,3,2]),[npln*nt,nd]);
tvis=kron(ones(nt,1),tlocal)+kron(npln*(0:nt-1)',0*tlocal+1);

if nargout<1
    %figure(1); clf;
    %set(axes,'FontSize',16);

    patch('vertices',pvis,'faces',tvis,'cdata',u(:), ...
               'facecol','interp','edgec','none');

    if nargin>=6 && ~isempty(clim)
      set(gca,'clim',clim);
    end
    %colorbar('FontSize',15,'Location','East');
    axis equal;axis tight;view(3);
end


% function surfaceplot(mesh,nref,UDG,ib,clim)
% 
% porder = mesh.porder;
% perm = mesh.perm;
% ne = mesh.ne;
% nf = mesh.nf;
% nd = mesh.nd;
% [npf,nfe] = size(perm);
% elcon = reshape(mesh.elcon,[npf nfe ne]);
% 
% in = ismember(mesh.f(:,end),ib);
% ns = sum(in); 
% dgnodes = zeros(npf,nd,ns);
% 
% sz = size(UDG);
% if length(sz)==2
%     u = reshape(UDG(1,:,:),[npf,nf]);
%     u = u(:,in);
% else
%     u = zeros(npf,ns);
% end
% 
% 
% % shapfc = mkshape(porder,mesh.plocfc,mesh.plocfc,mesh.elemtype);
% % dshapft = reshape(permute(shapfc(:,:,2:end),[2 3 1]),[npf*(nd-1) npf]);
% j = 1;
% for i = 1:nf
%     fi = mesh.f(i,end-1:end); % obtain two elements sharing the same face i  
% %     if fi(2)>0
% %         kf = mesh.t2f(fi(1),:);    % obtain neighboring faces 
% %         i1 = find(kf(1,:)==i);  % obtain the index of face i in the 1st element        
% %         j1 = elcon(:,i1,fi(1)) - (i-1)*npf;                
% %         kf = mesh.t2f(fi(2),:);    % obtain neighboring faces 
% %         i2 = find(kf(1,:)==i);  % obtain the index of face i in the second element        
% %         j2 = elcon(:,i2,fi(2)) - (i-1)*npf;                
% %         [mesh.dgnodes(perm(j1,i1),1:nd,fi(1))-mesh.dgnodes(perm(j2,i2),1:nd,fi(2))]                       
% %         pause
% %     end
%     if ismember(fi(2),ib) % face i is a boundary face
%         kf = mesh.t2f(fi(1),:);    % obtain neighboring faces 
%         i1 = find(kf(1,:)==i);  % obtain the index of face i in the 1st element        
%         j1 = elcon(:,i1,fi(1)) - (i-1)*npf;            
%         p = mesh.dgnodes(perm(j1,i1),1:nd,fi(1));        
%         dgnodes(:,:,j) = p;
%         u(:,j) = UDG(perm(j1,i1),1,fi(1));                
%         j = j + 1;                
% %         dpg = permute(reshape(dshapft*p,[npf nd-1 nd]),[1 3 2]);     
% %         nx = reshape(dpg(:,2,1,:).*dpg(:,3,2,:) - dpg(:,3,1,:).*dpg(:,2,2,:),[npf 1]);
% %         ny = reshape(dpg(:,3,1,:).*dpg(:,1,2,:) - dpg(:,1,1,:).*dpg(:,3,2,:),[npf 1]);
% %         nz = reshape(dpg(:,1,1,:).*dpg(:,2,2,:) - dpg(:,2,1,:).*dpg(:,1,2,:),[npf 1]);            
% %         detf = sqrt(nx.^2+ny.^2+nz.^2);
% %         nx = nx./detf;
% %         ny = ny./detf;
% %         nz = nz./detf;                        
% %         N = [nx ny nz];
% %         U = UDG(perm(j1,i1),:,fi(1));             
%     end
% end
% 
% plocal=mesh.plocfc;
% tlocal=mesh.tlocfc;
% if nref>0
%   if size(tlocal,2)==3  
%     A0=koornwinder(plocal(:,1:2),porder);
%     [plocal,tlocal]=uniref(plocal,tlocal,nref);
%     A=koornwinder(plocal(:,1:2),porder)/A0;
%   else
%     A0=tensorproduct(plocal(:,1:2),porder);
%     m = porder*(nref+1)+1;
%     [plocal,tlocal]=squaremesh(m,m,0,1);
%     A=tensorproduct(plocal(:,1:2),porder)/A0;  
%   end
%   npln=size(plocal,1);
%   sz=size(dgnodes); if length(sz)==2, sz = [sz,1]; end
%   dgnodes=reshape(A*reshape(dgnodes,npf,sz(2)*sz(3)),[npln,sz(2),sz(3)]);
%   u=A*u;
% end
% 
% nt = size(u,2);
% npln=size(plocal,1);
% nodesvis=reshape(permute(dgnodes,[1,3,2]),[npln*nt,nd]);
% tvis=kron(ones(nt,1),tlocal)+kron(npln*(0:nt-1)',0*tlocal+1);
% 
% figure(1); clf;
% set(axes,'FontSize',16);
% 
% patch('vertices',nodesvis,'faces',tvis,'cdata',u(:), ...
%            'facecol','interp','edgec','none');
%       
% if nargin>=5 && ~isempty(clim)
%   set(gca,'clim',clim);
% end
% colorbar('FontSize',15,'Location','East');
% axis equal;axis tight;view(3);
% 
% function e=boundedges(p,t,elemtype)
% %BOUNDEDGES Find boundary edges from triangular mesh
% %   E=BOUNDEDGES(P,T)
% 
% % Form all edges, non-duplicates are boundary edges
% 
% if elemtype==0
%     edges=[t(:,[1,2]);
%            t(:,[1,3]);
%            t(:,[2,3])];
%     node3=[t(:,3);t(:,2);t(:,1)];
% else
%     edges=[t(:,[1,2]);
%            t(:,[2,3]);
%            t(:,[3,4]);
%            t(:,[4,1]);];
%     node3=[t(:,4);t(:,3);t(:,2);t(:,1)];    
% end
% edges=sort(edges,2);
% [foo,ix,jx]=unique(edges,'rows');
% vec=histc(jx,1:max(jx));
% qx=find(vec==1);
% e=edges(ix(qx),:);
% node3=node3(ix(qx));
% 
% % Orientation
% v1=p(e(:,2),:)-p(e(:,1),:);
% v2=p(node3,:)-p(e(:,1),:);
% ix=find(v1(:,1).*v2(:,2)-v1(:,2).*v2(:,1)>0);
% e(ix,[1,2])=e(ix,[2,1]);
% 
% 
% function e1=segcollect(e)
% %SEGCOLLECT Collect polygons from edge segments.
% 
% ue=unique(e(:));
% he=histc(e(:),ue);
% current=ue(min(find(he==1))); % Find an endpoint
% if isempty(current) % Closed curve
%   current=e(1,1);
% end
% e1=current;
% while ~isempty(e)
%   ix=min(find(e(:,1)==e1(end)));
%   if isempty(ix)
%     ix=min(find(e(:,2)==e1(end)));
%     if isempty(ix) % >1 disjoint curves, recur
%       rest=segcollect(e);
%       e1={e1,rest{:}};
%       return;
%     end
%     next=e(ix,1);
%   else
%     next=e(ix,2);
%   end
%   e1=[e1,next];
%   e(ix,:)=[];
% end
% e1={e1};
% 
% 
% 
% 
