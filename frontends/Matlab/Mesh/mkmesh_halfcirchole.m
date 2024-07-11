function mesh=mkmesh_halfcirchole(porder)

R = 1;
W = 5;
L = 5;
H = 5;
n1= 40;
h2= 0.35;
growth=1.3;

phi=linspace(3*pi/2,pi/2,n1+1)';
phi=logdec(phi,1);
pv1=[R*cos(phi),R*sin(phi); 0,L; -W,L; -W,-H; 0 -H];

[mesh.p,mesh.t]=polymesh({pv1},1,[1,0],[h2,growth],@hh,[]);
[mesh.p,mesh.t] = fixmesh(mesh.p,mesh.t);
[mesh.f,mesh.t2f] = mkt2f(mesh.t);

bndexpr={sprintf('all(sqrt(sum(p.^2,2))<%g+1e-4)',R), ...         
         sprintf('all(p(:,2)>%g-1e-4)',L), ...
         sprintf('all(p(:,2)<-%g+1e-4)',H), ...
         sprintf('all(p(:,1)>-%g-1e-4)',0), ...          
         sprintf('all(p(:,1)<-%g+1e-4)',W)};
     
mesh.f = setbndnbrs(mesh.p,mesh.f,bndexpr);
mesh.fcurved = (mesh.f(:,4)==-1);
ic = mesh.fcurved;
mesh.tcurved = false(size(mesh.t,1),1);
mesh.tcurved(mesh.f(ic,3)) = true;

fd=@(p) sqrt(sum(p.^2,2))-R;
mesh.porder = porder;
%[mesh.plocal,mesh.tlocal] = uniformlocalpnts(mesh.porder);
%mesh.dgnodes = createnodes(mesh,fd);
[mesh.plocal,mesh.tlocal,mesh.plocfc,mesh.tlocfc,permnode,permedge,permface] = masternodes(mesh.porder,2,0,0);
mesh.permnode = permnode;
mesh.permedge = permedge;
mesh.permface = permface;
mesh.perm = permedge;

mesh.dgnodes = createnodes(mesh,fd);
mesh.ne = size(mesh.t,1);
mesh.np = size(mesh.p,1);
mesh.nf = size(mesh.f,1);
mesh.nd = 2;
mesh.elemtype = 0;
mesh.nodetype = 0;


function h=hh(p,hparabola)

h=inf+0*p(:,1);
if ~isempty(hparabola)
  p(:,1)=p(:,1)-hparabola(1);
  p=p(:,[2,1]);
  h1=hparabola(3)+hparabola(4)*abs(dparabola(p,hparabola(2)));
  h=min(h,h1);
end
