function boundaryplot(mesh,ib,color)
%MESHPLOT Plot mesh structure
%  MESHPLOT(MESH,[OPTS])
%
%  MESH:    Mesh structure
%  OPTS:    (logical)
%   OPTS(1): Plot elements/faces using p (default) or dgnodes
%   OPTS(2): Plot dgnodes
%   OPTS(3): Plot element numbers (2D only)
%   OPTS(4): Plot node numbers (2D only)
%   OPTS(5): Plot face numbers (2D only)
% if nargin<2 || isempty(opts), opts=0; end
% if length(opts)<2, opts=[opts,0]; end
% if length(opts)<3, opts=[opts,0]; end
% if length(opts)<4, opts=[opts,0]; end
% if length(opts)<5, opts=[opts,0]; end
%p=mesh.p';
%t=mesh.t';
%f=mesh.f';
dim=size(mesh.p,1);
%dpl=size(mesh.plocal,2);
% surface mesh
% if dpl==2 && dim==3
%   surface = 1;
% else
%   surface = 0;
% end
if dim < 1 || dim > 3
  error('Only can handle dim=1, dim=2 or dim=3');
end
pars={'facecolor',color,'edgecolor',color,'Linew',1,'FaceAlpha',1,'EdgeAlpha',1};
% if exist('hh','var')==0
%   hh=[];
% end
[nfe,~] = size(mesh.f);
[indf, inde] = find(mesh.f==ib);
nve = size(mesh.t,1);
nvf = dim;
if (dim==3) && (nfe==6)
 nvf = 4;
end
elemtype = 0;
if (dim==2) && (nve==4)
 elemtype = 1;
elseif (dim==3) && (nve==8)
 elemtype = 1;
end
face = getelemface(dim,elemtype);
nf = length(inde);
bf = zeros(nf, nvf);
for i = 1:nf
 ei = inde(i);
 fi = indf(i);
 bf(i,:) = mesh.t(face(:,fi),ei);
end
% figure(1); clf;
patch('faces',bf,'vertices',mesh.p',pars{:});
% if opts(1)==0 % plot boundary faces using p
%   hh=[hh;patch('faces',bf,'vertices',p,pars{:})];
% else     % plot boundary faces using dgnodes
%   %bf = 1:size(f,1);
%   e=boundedges(mesh.plocfc,mesh.tlocfc,mesh.elemtype);
%   e1=segcollect(e);
%   axis equal,axis off
%   nt=length(bf);
%   hh=zeros(nt,1);
%   for it=1:nt
%     el = bf(bf(it),end-1);
%     fc = mesh.t2f(el,:);
%     fi = find(fc==bf(it));
%     px=mesh.dgnodes(mesh.perm(:,fi),1,el);
%     py=mesh.dgnodes(mesh.perm(:,fi),2,el);
%     pz=mesh.dgnodes(mesh.perm(:,fi),3,el);
%     %if mean(pz(:))>=45
%     %pw=udg(mesh.perm(:,fi),1,el);
%     %if max(abs(py))<1e-10
%     %hh(it)=patch(px(e1{1}'),py(e1{1}'),pz(e1{1}'),pw(e1{1}'),pars{:});
%     hh(it)=patch(px(e1{1}'),py(e1{1}'),pz(e1{1}'),0*pz(e1{1}'),pars{:});
%     %hh(it)=patch(px(e1{1}'),py(e1{1}'),pz(e1{1}'),pw(e1{1}'));
%     %end
%     %hh(it)=patch(px,py,pz,0.0*px,pars{:});
%     %end
%   end
% end
if dim==3
view(3),axis equal; % colorbar;
end
% function e=boundedges(p,t,elemtype)
% %BOUNDEDGES Find boundary edges from triangular mesh
% %  E=BOUNDEDGES(P,T)
%
% % Form all edges, non-duplicates are boundary edges
%
% if elemtype==0
%   edges=[t(:,[1,2]);
%      t(:,[1,3]);
%      t(:,[2,3])];
%   node3=[t(:,3);t(:,2);t(:,1)];
% else
%   edges=[t(:,[1,2]);
%      t(:,[2,3]);
%      t(:,[3,4]);
%      t(:,[4,1]);];
%   node3=[t(:,4);t(:,3);t(:,2);t(:,1)];
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
%  current=e(1,1);
% end
% e1=current;
% while ~isempty(e)
%  ix=min(find(e(:,1)==e1(end)));
%  if isempty(ix)
%   ix=min(find(e(:,2)==e1(end)));
%   if isempty(ix) % >1 disjoint curves, recur
%    rest=segcollect(e);
%    e1={e1,rest{:}};
%    return;
%   end
%   next=e(ix,1);
%  else
%   next=e(ix,2);
%  end
%  e1=[e1,next];
%  e(ix,:)=[];
% end
% e1={e1};
%
%
%