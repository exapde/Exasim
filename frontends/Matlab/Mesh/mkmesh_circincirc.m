function mesh = mkmesh_circincirc(porder,m,n,r1,r2)
% Exampe: mesh = mkmesh_circincirc(2,21,21,1,10);

nodetype = 1;
elemtype = 1;
[p,t] = squaremesh(m,n,1,elemtype);
p(:,2) = loginc(p(:,2),0.2);

bndexpr = {'all(p(:,2)<1e-3)','all(p(:,1)>max(p0(:,1))-1e-3)', ...
           'all(p(:,2)>max(p0(:,2))-1e-3)','all(p(:,1)<1e-3)'};     
mesh = mkmesh(p,t,porder,bndexpr,elemtype,nodetype);

p = mesh.p;
pnew(:,1) = (r1+(r2-r1)*p(:,2)).*sin(2*pi*p(:,1));
pnew(:,2) = (r1+(r2-r1)*p(:,2)).*cos(2*pi*p(:,1));
mesh.p = pnew;
[mesh.p,mesh.t]=fixmesh(mesh.p,mesh.t);

p = mesh.dgnodes;   
pnew = zeros(size(p));
pnew(:,1,:) = (r1+(r2-r1)*p(:,2,:)).*sin(2*pi*p(:,1,:));
pnew(:,2,:) = (r1+(r2-r1)*p(:,2,:)).*cos(2*pi*p(:,1,:));
mesh.dgnodes = pnew;
mesh.fcurved = true(size(mesh.f,1),1);
mesh.tcurved = true(size(mesh.t,1),1);

bndexpr = {'all(sqrt(sum(p.^2,2))<1.15)','true'};     
[mesh.f,mesh.t2f] = mkt2f(mesh.t,elemtype);
mesh.f = setbndnbrs(mesh.p,mesh.f,bndexpr);
mesh.nf = size(mesh.f,1);
