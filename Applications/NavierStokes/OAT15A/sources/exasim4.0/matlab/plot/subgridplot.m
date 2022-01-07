function hh = subgridplot(mesh,mastersubgrid,opts)
%DGMESHPLOT plot the multiscale DG mesh
%
%   DGMESHPLOT(mesh,mastersubgrid,opts)
%
%    MESH:       Mesh structure
%    MASTERGRID: Mastergrid structure
%    OPTS:       (logical)
%      OPTS(1):  Plot dgnodes 
%      OPTS(2):  Plot element numbers (2D only)
%      OPTS(3):  Plot node numbers (2D only)
%      OPTS(4):  Plot face numbers (2D only)
%

if nargin<3 || isempty(opts), opts=0; end
if length(opts)<2, opts=[opts,0]; end
if length(opts)<3, opts=[opts,0]; end
if length(opts)<4, opts=[opts,0]; end

% dimension
nd = mesh.nd;
dim= nd;

% number of subgrids
ns = max(mesh.subgrids);

dpl=size(mesh.plocal,2);

% surface mesh
if dpl==2 && dim==3
    surface = 1;
else 
    surface = 0;
end
    
if dim ~= 2
    error('Only can handle dim = 2');
end

pars={'facecolor',[0.2,0.5,1.0],'edgecolor','r','Linew',0.5};

if exist('hh','var')==0
    hh=[];
end

for j = 1:ns % for each subgrid
    % obtain element indices for subgrid j
    ind = find(mesh.subgrids==j);    
    nj  = length(ind);
    
    % obtain geometry nodes
    dgnodes = mesh.dgnodes(:,:,ind);                
    
    if dim == 2 || surface==1    
        e=boundedges(mesh.plocal,mesh.tlocal,mesh.elemtype);
        e1=segcollect(e);        
        axis equal,axis off;

        % obtain subgrid geometry nodes
        geomnodes = mastersubgrid{j}.geomnodes;        
        ng = size(geomnodes,3);

        % compute shape functions at the subgrid DG nodes
        tmp = mkshape(mesh.porder,mesh.plocal,geomnodes(:,:,1),mesh.elemtype);
        n1 = size(tmp,1); 
        n2 = size(tmp,2); 
        shp = zeros(n2,n1,ng);                
        for k=1:ng
            tmp = mkshape(mesh.porder,mesh.plocal,geomnodes(:,:,k),mesh.elemtype);
            shp(:,:,k) = tmp(:,:,1)';
        end      

        % compute subgrid DG nodes in the physical space
        dgx = zeros(n2,nd,nj*ng);
        for i=1:nj % for each superelement                           
            for k=1:ng % for each subgrid element                    
                dgx(:,:,(i-1)*ng+k) = shp(:,:,k)*dgnodes(:,:,i);                   
            end
        end        

        nt=size(dgx,3);
        hh=zeros(nt,1);
        for it=1:nt
            px=dgx(:,1,it);
            py=dgx(:,2,it);
            if surface==1  
                pz=dgx(:,3,it);
            else
                pz=0*px;
            end
            hh(it)=patch(px(e1{1}'),py(e1{1}'),pz(e1{1}'),0.0*e1{1}',pars{:});
        end                
        if surface==0
            view(2),axis equal;
        else
            view(3),axis equal;
        end            
    end    
    hold on;        
    
    if opts(1)==1
        % obtain subgrid DG nodes
        elemnodes = mastersubgrid{j}.elemnodes;  
        ng = size(elemnodes,3);
        tmp = mkshape(mesh.porder,mesh.plocal,elemnodes(:,:,1),mesh.elemtype);
        n1 = size(tmp,1); 
        n2 = size(tmp,2); 
        shp = zeros(n2,n1,ng);                            
        for k=1:ng
            tmp = mkshape(mesh.porder,mesh.plocal,elemnodes(:,:,k),mesh.elemtype);
            shp(:,:,k) = tmp(:,:,1)';
        end        

        % compute subgrid DG nodes in the physical space
        dgx = zeros(n2,nd,nj*ng);
        for i=1:nj % for each element                 
            for k=1:ng % for each subgrid element
                dgx(:,:,(i-1)*ng+k) = shp(:,:,k)*dgnodes(:,:,i);                    
            end
        end                 
        xx=squeeze(dgx(:,1,:));
        yy=squeeze(dgx(:,2,:));            
        line(xx(:),yy(:),0*xx(:),'lines','n','marker','.','markersize',16,'col','b');    
    end            
end

p = mesh.p;
t = mesh.t;

if opts(2)
    if dim == 2 || surface==1   
        for it=1:size(t,1)
            pmid=mean(p(t(it,:),:),1);
            txtpars={'fontname','times','fontsize',20,'horizontala','center','BackgroundColor',[1,1,1]};
            text(pmid(1),pmid(2),num2str(it),txtpars{:});
        end
    end
end

if opts(3)
    if dim == 2 || surface==1
        for i=1:size(p,1)
            txtpars={'fontname','times','fontsize',20,'fontweight','bold', ...
                'horizontala','center','col','w', 'BackgroundColor',[0.5,0.5,0.5]};
            text(p(i,1),p(i,2),num2str(i),txtpars{:});
        end
    end
end
        
if nargout<1, clear hh; end

set(gcf,'rend','z');
axis equal;
drawnow

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





