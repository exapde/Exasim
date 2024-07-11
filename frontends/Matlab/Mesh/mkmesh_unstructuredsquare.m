function mesh = mkmesh_unstructuredsquare(h,porder)
%MKMESH_SQUARE Creates 2D tri/quad mesh data structure for unit square.
%   MESH=MKMESH_SQUARE(M,N,PORDER,PARITY)
%
%      MESH:      Mesh structure
%      M:         Number of points in the horizaontal direction 
%      N:         Number of points in the vertical direction
%      PORDER:    Polynomial Order of Approximation (default=1)
%      PARITY:    Flag determining the the triangular pattern
%                 Flag = 0 (diagonals SW - NE) (default)
%                 Flag = 1 (diagonals NW - SE)
%      ELEMTYPE:  Flag determining element type
%                 Flag = 0 tri/tet elements (default)
%                 Flag = 1 quad/hex elements
%      NODETYPE:  Flag determining node distribution 
%                 Flag = 0 uniform distribution (default)
%                 Flag = 1 nonuniform distribution

%   See also: SQUAREMESH, MKMESH
%

if nargin<1, h=1/8; end
if nargin<2, porder=1; end

%L = 1.5;
pv1 = [0 0; 0 0.75; 1 0.75; 1 0];
[p1,t1]=polymesh({pv1},[1],[0,1],[h/1.2,1.3]);
ind = find(p1(:,2)==0.75);
ind = [ind(1); ind(3:end); ind(2)];
pv2 = [p1(ind,1) p1(ind,2); 1 1; 0 1];
[p2,t2]=polymesh({pv2},[1],[1,0],[h/1.2,1.3]);

[p,t] = connectmesh(p1,t1,p2,t2);

%figure(1); clf; simpplot(p,t); axis on; axis tight;

bndexpr = {'all(p(:,2)<1e-3)','all(p(:,1)>max(p0(:,1))-1e-3)', ...
           'all(p(:,2)>max(p0(:,2))-1e-3)','all(p(:,1)<1e-3)'};     

mesh = mkmesh(p,t,porder,bndexpr,0,1);

pmid=(p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;
ne = size(pmid,1);
subgrids = ones(ne,1);
for n=1:ne
    if pmid(n,2)>0.48 && abs(pmid(n,1)-0.5)<0.05;
        subgrids(n) = 2;
    end
%     if pmid(n,2)>0.75 && abs(pmid(n,1)-0.5)<0.06;
%         subgrids(n) = 2;
%     end
end
mesh.subgrids = subgrids;


ind1 = subgrids==1;
ind2 = subgrids==2;

bcol=[.8,1,.8];
figure(1); clf;
trimesh(t(ind1,:),p(:,1),p(:,2),0*p(:,1),'facecolor',bcol,'edgecolor','k');
hold on;
trimesh(t(ind2,:),p(:,1),p(:,2),0*p(:,1),'facecolor',[1 0 0],'edgecolor','k');

view(2)
axis equal
axis on
ax=axis;axis(ax*1.001);

% pmid=(p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;
% ne = size(pmid,1);
% subgrids = ones(ne,1);
% for n=1:ne
%     if pmid(n,2)>0.485 && abs(pmid(n,1)-0.5)<0.04;
%         subgrids(n) = 3;
%     elseif pmid(n,2)>0.475 && abs(pmid(n,1)-0.5)<0.08;
%         subgrids(n) = 2;
%     end
% end
% mesh.subgrids = subgrids;
% 
% 
% ind1 = subgrids==1;
% ind2 = subgrids==2;
% ind3 = subgrids==3;
% 
% bcol=[.8,1,.8];
% figure(1); clf;
% trimesh(t(ind1,:),p(:,1),p(:,2),0*p(:,1),'facecolor',bcol,'edgecolor','k');
% hold on;
% trimesh(t(ind2,:),p(:,1),p(:,2),0*p(:,1),'facecolor',[0 0 1],'edgecolor','k');
% trimesh(t(ind3,:),p(:,1),p(:,2),0*p(:,1),'facecolor',[1 0 0],'edgecolor','k');
% 
% view(2)
% axis equal
% axis on
% ax=axis;axis(ax*1.001);
