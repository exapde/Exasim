function ssplot(mesh,mastersubgrid,u,iu,clim,nref)

% number of subgrids
ns = max(mesh.subgrids);
nd = mesh.nd;
mesht = mesh;


cla;
for j = 1:ns % for each subgrid
    % obtain element indices for subgrid j
    ind = find(mesh.subgrids==j);    
    if isempty(ind)==0
        nj = length(ind);
%         ng = size(mastersubgrid{j}.elemnodes,3);
%         npm = size(mesh.dgnodes,1);        
%         npv = size(u{j},1);
        dgnodes = mesh.dgnodes(:,:,ind);                                        
%         udg = subgridprojection(mastersubgrid{j},mastersubgrid{1},mesht.dgnodes,reshape(u{j}(:,iu,:),[npv 1 ng nj]));                
%         udg = reshape(udg,[npm 1 nj]);        
%         scalarplot(mesht,udg,nref);        
        
        % obtain subgrid interpolation nodes
        if mastersubgrid{j}.porder==0
            [plocvl,tlocvl] = masternodes(1,nd,mastersubgrid{j}.elemtype,mastersubgrid{j}.nodetype);
            elemnodes = mknodes(mastersubgrid{j}.p,mastersubgrid{j}.t,plocvl);            
            mesht.plocal  = plocvl;
            mesht.tlocal  = tlocvl;                
        else
            elemnodes = mastersubgrid{j}.elemnodes;    
            mesht.plocal  = mastersubgrid{j}.plocvl;
            mesht.tlocal  = mastersubgrid{j}.tlocvl;                
        end
%         if mastersubgrid{j}.porder==0
%             qorder = 1;
%             [plocvl,tlocvl] = masternodes(qorder,nd,mastersubgrid{j}.elemtype,mastersubgrid{j}.nodetype);
%             elemnodes = mknodes(mastersubgrid{j}.p,mastersubgrid{j}.t,plocvl);  
%             mesht.plocal  = plocvl;
%             mesht.tlocal  = tlocvl;                
%         else
%             elemnodes = mastersubgrid{j}.geomnodes;
%         end
        ng = size(elemnodes,3);

        % compute shape functions at the subgrid DG nodes
%         tmp = mkshape(mesh.porder,mesh.plocal,elemnodes(:,:,1),mesh.elemtype);
%         n1 = size(tmp,1); 
%         n2 = size(tmp,2); 
%         shp = zeros(npm,n1,ng);                
%         shv = zeros(n2,n0,ng);
        clear shp shv;
        for k=1:ng
            tmp = mkshape(mesh.porder,mesh.plocal,elemnodes(:,:,k),mesh.elemtype);
            shp(:,:,k) = tmp(:,:,1)';
            tmp = mkshape(mastersubgrid{j}.porder,mastersubgrid{j}.plocvl,elemnodes(:,:,k),mesh.elemtype);
            shv(:,:,k) = tmp(:,:,1)';
        end            
                
        % compute subgrid DG nodes in the physical space
        %nc  = size(u{j},2);        
        npm = size(shp,1);        
        dgx = zeros(npm,nd,nj*ng);        
        udg = zeros(npm,1,nj*ng);        
        for i=1:nj % for each superelement                           
            for k=1:ng % for each subgrid element                    
                dgx(:,:,(i-1)*ng+k) = shp(:,:,k)*dgnodes(:,:,i);                   
                udg(:,1,(i-1)*ng+k) = shv(:,:,k)*u{j}(:,iu,(i-1)*ng+k);                   
            end
        end        
        
        mesht.dgnodes = dgx; 
        mesht.t       = mesh.t(ind,:);
        %mesht.porder  = mastersubgrid{j}.porder;   
        %npv = size(u,1);
        %nc  = size(u,2);
        %nes = size(mastersubgrid{j}.t,1);        
        %scalarplot(mesht,u{j}(:,:,(ind(1)-1)*nes+1:ind(end)*nes),nref)        
%         size(dgx)
%         size(u{j})
        scalarplot(mesht,udg,nref);
        hold on;
    end
end

if nargin>=3 && ~isempty(clim)
  set(gca,'clim',clim);
end

trimesh(mesh.t,mesh.p(:,1),mesh.p(:,2),0*mesh.p(:,1),'facecolor','none','edgecolor','k');
set(gcf,'rend','z');
colorbar;
colorbar('FontSize',15);
axis equal;
axis tight;
box on;
axis off;
drawnow

% plot the mesh


% if nargin>=6 && ~isempty(surf)
%    cameramenu;
% end

function scalarplot(mesh,u,nref)
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
porder=mesh.porder;

% if porder==0
%   nref=0; 
%   u=reshape(u,1,nt);
%   u=repmat(u,[npl,1]);
% else
%   u=reshape(u,npl,nt);
% end
u=reshape(u,npl,nt);
plocal=mesh.plocal;
tlocal=mesh.tlocal;  
dgnodes=mesh.dgnodes;
dgnodes = dgnodes(:,1:nd,:);

if nargin<3 || isempty(nref), nref=ceil(log2(max(porder,1))); end

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

patch('vertices',nodesvis,'faces',tvis,'cdata',u, ...
      'facecol','interp','edgec','none');
      
    % Plot curved mesh
% e=boundedges(nodesvis,tvis);
% dgnodesltx=nodesvis(:,1);
% dgnodeslty=nodesvis(:,2);
% line(dgnodesltx(e'),dgnodeslty(e'),'color',[0,0,0],'LineWidth',1);    

% 
% function e=boundedges(p,t)
% %BOUNDEDGES Find boundary edges from triangular mesh
% %   E=BOUNDEDGES(P,T)
% 
% % Form all edges, non-duplicates are boundary edges
% edges=[t(:,[1,2]);
%        t(:,[1,3]);
%        t(:,[2,3])];
% node3=[t(:,3);t(:,2);t(:,1)];
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
