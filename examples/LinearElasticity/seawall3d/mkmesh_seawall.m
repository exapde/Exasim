function mesh = mkmesh_seawall(porder, elemtype)

[pde,mesh] = initializeexasim();

[wall, land] = lbseawall();
poly2gmsh('seawall.geo', wall, elemtype, 0.5);
[mesh.p, mesh.t] = gmshcall(pde, "seawall", 2, elemtype);
mesh.t = mesh.t(end:-1:1,:);

figure(1);clf;simpplot(mesh.p',mesh.t');
hold on;
%plot(wall(:,1), wall(:,2), 'o');
axis on; axis equal;

bndexpr = {'true'}; 
mesh = mkmesh(mesh.p',mesh.t',porder,bndexpr,elemtype,1);
mesh.p = mesh.p';
mesh.t = mesh.t';

ymin = min(mesh.p(2,:));
ymax = max(mesh.p(2,:));

mesh.boundaryexpr = {@(p) abs(p(2,:)-ymin)<1e-3, @(p) abs(p(2,:)-ymax)<1e-3, @(p) p(1,:)<1e-3, @(p) abs(p(1,:)-land(1,1))<1e-3, @(p) abs(p(2,:))<100};
mesh.boundarycondition = [1;1;3;2]; 
mesh.periodicexpr = [];
mesh.f = facenumbering(mesh.p,mesh.t,elemtype,mesh.boundaryexpr,mesh.periodicexpr);

%plotboundary(mesh);