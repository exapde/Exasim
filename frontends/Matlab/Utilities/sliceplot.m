function [pp,tt,cc] = sliceplot(mesh, udg, c, nref)

if mesh.nd ~= 3
    error('isosurface plot is only for 3D field.');
end
if (mesh.porder == 0) || (size(udg,1) == 1)
    error('porder must be greater than 0.');
end

if nargin <= 4
    nref = [];   
end

nd = mesh.nd;
[ne,nv] = size(mesh.t);

p = reshape(mesh.p(mesh.t',:),[nv ne nd]);
d = surfunc(p,c);

% elements contain the surface
inde = ~(all(d>=0,1) | all(d<=0,1));
inde = find(inde==1);

if isempty(inde)
    warning('No surface exists');
    return;
end

npl = size(udg,1);
nel = length(inde);

% get scalar fields and coordinates on those elements
udg = udg(:,inde);
pdg = permute(mesh.dgnodes(:,1:nd,inde),[1 3 2]);
 
if mesh.porder>1
    % subdivision 
    [pdg,udg] = scalarrefine(mesh,pdg,udg,nref);    
    
    [npl,nel] = size(udg);    
end

% if hexes then convert them into tets
if npl==8
    % cube to tets
    c2t = [1 2 4 6; 1 5 4 6; 5 8 4 6; 2 3 4 6; 7 3 4 6; 7 8 4 6]';

    udg = reshape(udg(c2t,:),[4 6*nel]);    
    pdg = reshape(pdg(c2t,:),[4 6*nel nd]);    
    
    nel = 6*nel;
end

d1 = reshape(surfunc(pdg(1,:,:),c),[nel 1]);
d2 = reshape(surfunc(pdg(2,:,:),c),[nel 1]);
d3 = reshape(surfunc(pdg(3,:,:),c),[nel 1]);
d4 = reshape(surfunc(pdg(4,:,:),c),[nel 1]);
d1d2 = d1.*d2;
d1d3 = d1.*d3;
d1d4 = d1.*d4;
d2d3 = d2.*d3;
d2d4 = d2.*d4;
d3d4 = d3.*d4;

% case 1
in1 = find( sum(abs(d1)+abs(d2)+abs(d3))<1e-10 );
ne1 = length(in1);
pp{1} = pdg(1:nd,in1,:);
tt{1} = reshape(1:nd*ne1,[nd,ne1])';
cc{1} = udg(1:nd,in1);
pp{1} = reshape(pp{1},[nd*ne1,nd]);

% case 2
in2 = find( sum(abs(d1)+abs(d2)+abs(d4))<1e-10 );
ne2 = length(in2);
pp{2} = pdg([1 2 4],in2,:);
tt{2} = reshape(1:nd*ne2,[nd,ne2])';
cc{2} = udg([1 2 4],in2);
pp{2} = reshape(pp{2},[nd*ne2,nd]);

% case 3
in3 = find( sum(abs(d1)+abs(d3)+abs(d4))<1e-10 );
ne3 = length(in3);
pp{3} = pdg([1 3 4],in3,:);
tt{3} = reshape(1:nd*ne3,[nd,ne3])';
cc{3} = udg([1 3 4],in3);
pp{3} = reshape(pp{3},[nd*ne3,nd]);

% case 4
in4 = find( sum(abs(d2)+abs(d3)+abs(d4))<1e-10 );
ne4 = length(in4);
pp{4} = pdg([2 3 4],in4,:);
tt{4} = reshape(1:nd*ne4,[nd,ne4])';
cc{4} = udg([2 3 4],in4);
pp{4} = reshape(pp{4},[nd*ne4,nd]);

% case 5
in5 = find( d1d2<=0 & d1d3<=0 & d1d4 <=0 );
ne5 = length(in5);
pp{5} = zeros(nd,ne5,nd);
tt{5} = reshape(1:nd*ne5,[nd,ne5])';
cc{5} = zeros(nd,ne5);
for j = 1:nd
    k = j + 1;
    [p,a] = InterpolateVertices(pdg(1,in5,:),pdg(k,in5,:),c);
%     size(p)
%     size(pp{5})
    pp{5}(j,:,:) = p;
    cc{5}(j,:,:) = a.*udg(1,in5) + (1-a).*udg(k,in5);    
end
y = pp{5}(:,:,2);

% case 6
in6 = find( d1d2<=0 & d2d3<=0 & d2d4 <= 0 );
ne6 = length(in6);
pp{6} = zeros(nd,ne6,nd);
tt{6} = reshape(1:nd*ne6,[nd,ne6])';
cc{6} = zeros(nd,ne6);
for j = 1:nd
    if j == 1
        k = j;
    else
        k = j + 1;
    end
    [p,a] = InterpolateVertices(pdg(2,in6,:),pdg(k,in6,:),c);
    pp{6}(j,:,:) = p;
    cc{6}(j,:,:) = a.*udg(2,in6) + (1-a).*udg(k,in6);
end

% case 7
in7 = find( d1d3<=0 & d2d3<=0 & d3d4 <= 0 );
ne7 = length(in7);
pp{7} = zeros(nd,ne7,nd);
tt{7} = reshape(1:nd*ne7,[nd,ne7])';
cc{7} = zeros(nd,ne7);
for j = 1:nd
    if j < 3
        k = j;
    else
        k = j + 1;
    end
    [p,a] = InterpolateVertices(pdg(3,in7,:),pdg(k,in7,:),c);
    pp{7}(j,:,:) = p;
    cc{7}(j,:,:) = a.*udg(3,in7) + (1-a).*udg(k,in7);
end

% case 8
in8 = find( d1d4<=0 & d2d4<=0 & d3d4 <=0 );
ne8 = length(in8);
pp{8} = zeros(nd,ne8,nd);
tt{8} = reshape(1:nd*ne8,[nd,ne8])';
cc{8} = zeros(nd,ne8);
for j = 1:nd    
    k = j;    
    [p,a] = InterpolateVertices(pdg(4,in8,:),pdg(k,in8,:),c);
    pp{8}(j,:,:) = p;
    cc{8}(j,:,:) = a.*udg(4,in8) + (1-a).*udg(k,in8);
end

% case 9
in9 = find( d1d2<=0 & d1d4<=0 & d2d3 <= 0 & d3d4 <= 0 );
ne9 = length(in9);
nv = 4;
pp{9} = zeros(nv,ne9,nd);
tt{9} = reshape(1:nv*ne9,[nv,ne9])';
cc{9} = zeros(nv,ne9);
k1 = [1 2 3 1];
k2 = [2 3 4 4];
for j = 1:nv
    [p,a] = InterpolateVertices(pdg(k1(j),in9,:),pdg(k2(j),in9,:),c);
    pp{9}(j,:,:) = p;
    cc{9}(j,:,:) = a.*udg(k1(j),in9) + (1-a).*udg(k2(j),in9);
end

% case 10
in10 = find( d1d2<=0 & d1d3<=0 & d2d4 <= 0 & d3d4 <= 0 );
ne10 = length(in10);
pp{10} = zeros(nv,ne10,nd);
tt{10} = reshape(1:nv*ne10,[nv,ne10])';
cc{10} = zeros(nv,ne10);
k1 = [1 1 3 2];
k2 = [2 3 4 4];
for j = 1:nv
    [p,a] = InterpolateVertices(pdg(k1(j),in10,:),pdg(k2(j),in10,:),c);
    pp{10}(j,:,:) = p;
    cc{10}(j,:,:) = a.*udg(k1(j),in10) + (1-a).*udg(k2(j),in10);
end

% case 11
in11 = find( d1d3<=0 & d1d4<=0 & d2d3 <= 0 & d2d4 <= 0 );
ne11 = length(in11);
pp{11} = zeros(nv,ne11,nd);
tt{11} = reshape(1:nv*ne11,[nv,ne11])';
cc{11} = zeros(nv,ne11);
k1 = [1 1 2 2];
k2 = [3 4 4 3];
for j = 1:nv
    [p,a] = InterpolateVertices(pdg(k1(j),in11,:),pdg(k2(j),in11,:),c);
    pp{11}(j,:,:) = p;
    cc{11}(j,:,:) = a.*udg(k1(j),in11) + (1-a).*udg(k2(j),in11);
end

i = 5; pp{i} = reshape(pp{i},[nd*ne5,nd]);
i = 6; pp{i} = reshape(pp{i},[nd*ne6,nd]);
i = 7; pp{i} = reshape(pp{i},[nd*ne7,nd]);
i = 8; pp{i} = reshape(pp{i},[nd*ne8,nd]);
i = 9; pp{i} = reshape(pp{i},[nv*ne9,nd]);
i = 10; pp{i} = reshape(pp{i},[nv*ne10,nd]);
i = 11; pp{i} = reshape(pp{i},[nv*ne11,nd]);

%figure(1);clf;
hold on;
for i = 1:length(cc)
    if isempty(pp{i})==0
    patch('vertices',pp{i},'faces',tt{i},'cdata',cc{i}(:), ...
       'facecol','interp','edgec','none','FaceLighting','gouraud'); 
    end
end
hold off;    
colormap jet;    
lighting gouraud
xlabel('x');
ylabel('y');
zlabel('z');
xmin = min(mesh.p(:,1)); xmax = max(mesh.p(:,1));
ymin = min(mesh.p(:,2)); ymax = max(mesh.p(:,2));
zmin = min(mesh.p(:,3)); zmax = max(mesh.p(:,3));
axis([xmin xmax ymin ymax zmin zmax]);
[xmin xmax ymin ymax zmin zmax]
colorbar('FontSize',12);
set(gca,'FontSize',16);
%set(gca,'clim',[-0.7 0.7]);
% cmin = zeros(1,length(cc)); cmax = cmin;
% for i = 1:length(cc)
%     cmin(i) = min(cc{i}(:)); cmax(i) = max(cc{i}(:));
% end
% set(gca,'clim',[min(cmin) max(cmax)]);

% cameratoolbar('SetCoordSys','y');
% campos([-2 5 5]);
% camtarget([0.65 0.325 0.05]);
% camzoom(1.5);
% camlight('right');    
% axis off;

function [pref,uref] = scalarrefine(mesh,p,u,nref)

[npl, nt, nd] = size(p);
porder=mesh.porder;
plocal=mesh.plocal;
tlocal=mesh.tlocal;

if isempty(nref), nref=ceil(log2(max(porder,1))); end
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
t = kron(ones(nt,1),tlocal)+kron(npln*(0:nt-1)',0*tlocal+1);
ne = size(t,1);
np = npln*nt;

uref = reshape(A*u,[np 1]);
uref = reshape(uref(t'), [size(t,2) ne]);

pref = reshape(A*reshape(p,npl,nt*nd),[np,nd]);
pref = reshape(pref(t',:),[size(t,2) ne nd]);


function d = surfunc(p,c)

d = c(end);
s = ndims(p);
if s==1
    for i = 1:length(c)-1
        d = d + c(i)*p(i);
    end
elseif s==2
    for i = 1:length(c)-1
        d = d + c(i)*p(:,i);
    end   
elseif s==3
    for i = 1:length(c)-1
        d = d + c(i)*p(:,:,i);
    end
end


function [p,a] = InterpolateVertices(p1,p2,c)

c0 = c; c0(end) = 0;
a = -surfunc(p2,c)./(surfunc(p1,c0)-surfunc(p2,c0));

p = 0*p1;
for i = 1:length(c)-1
    p(:,:,i) = a.*p1(:,:,i) + (1-a).*p2(:,:,i);
end


