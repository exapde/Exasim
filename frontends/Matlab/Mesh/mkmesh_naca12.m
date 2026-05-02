function mesh = mkmesh_naca12(porder)

[p,t] = naca0012();

elemtype = 0;
nodetype = 1;
bndexpr = {'all(sqrt(sum(p.^2,2))<2)','all(sqrt(sum(p.^2,2))>2)'};  

mesh = mkmesh(p,t,porder,bndexpr,elemtype,nodetype);

fb = 1;
mesh.fcurved = (mesh.f(:,end)==-fb);
ic = mesh.fcurved;
mesh.tcurved = false(size(mesh.t,1),1);
mesh.tcurved(mesh.f(ic,end-1)) = true;
mesh.dgnodes = makedgnodes(mesh,@naca12dist,fb);

mesh.p = mesh.p';
mesh.t = mesh.t';
mesh.boundaryexpr = {@(p) sqrt((p(1,:)-.5).^2+p(2,:).^2)<3, @(p) abs(p(1,:))<20};
mesh.boundarycondition = [1;2];
mesh.curvedboundary = [1 0];
mesh.curvedboundaryexpr = {@(p) p(2,:).^2-(5*0.01*12*(0.29690*sqrt(abs(p(1,:)))-0.12600*p(1,:)-0.35160*p(1,:).^2+0.28430*p(1,:).^3-0.10150*p(1,:).^4)).^2, @(p) 0};
mesh.periodicexpr =  {};
mesh.f = facenumbering(mesh.p,mesh.t,mesh.elemtype,mesh.boundaryexpr,mesh.periodicexpr);
mesh.xpe = mesh.plocal; 
mesh.telem = mesh.tlocal;

figure(2);clf;boundaryplot(mesh,1); hold on; %inlet
boundaryplot(mesh,2); %wall

