function mesh=mkmesh_qcirchole(porder,h1,h2,hparabola,width,height,growth)

if nargin<1, porder=1; end
if nargin<2, h1=0.4; end
if nargin<3, h2=1.5; end
if nargin<4, hparabola=[]; end % [-1.55,.2,.02,.5]
if nargin<6, width=8; end
if nargin<6, height=10; end
if nargin<7, growth=1.3; end

n1=ceil(pi/2/h1);
phi=linspace(pi,pi/2,n1+1)';
phi=logdec(phi,1);
pv1=[cos(phi),sin(phi); 0,height; -width,height; -width,0];

[mesh.p,mesh.t]=polymesh({pv1},[1],[1,0],[h2,growth],@hh,hparabola);
[mesh.p,mesh.t] = fixmesh(mesh.p,mesh.t);
[mesh.f,mesh.t2f] = mkt2f(mesh.t);

bndexpr={'all(p(:,2)<1e-3)',...
         'all(sqrt(sum(p.^2,2))<1+1e-3)', ...
         'all(p(:,1)>-1e-3)', ...
         sprintf('all(p(:,2)>%g-1e-3)',height), ...
         sprintf('all(p(:,1)<-%g+1e-3)',width)};
     
mesh.f = setbndnbrs(mesh.p,mesh.f,bndexpr);
mesh.fcurved = (mesh.f(:,4)==-2);
ic = mesh.fcurved;
mesh.tcurved = false(size(mesh.t,1),1);
mesh.tcurved(mesh.f(ic,3)) = true;

fd=@(p) sqrt(sum(p.^2,2))-1;
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
