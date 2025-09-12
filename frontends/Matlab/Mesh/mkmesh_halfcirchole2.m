function mesh=mkmesh_halfcirchole2(porder)

if nargin<2, mul1=1; end
if nargin<3, mul2=1; end

xmin = -10;
xmax = 10;
ymin = 0;
ymax = 10;

m=4*mul1;
n=12*mul1;
o=4*mul2;

s1=(1:m)/m; z1=0*s1;
s2=(1:n)/n; z2=0*s2;
%phi=2*pi*(1:4*o)/(4*o);
phi=linspace(0,pi,4*o);

% pv1=[xmin+(xmax-xmin)*s2,xmax+z1,xmax-(xmax-xmin)*s2,xmin+z1;
%      ymin+z2,ymin+(ymax-ymin)*s1,ymax+z2,ymax-(ymax-ymin)*s1]';
pv1=[xmin xmin xmax xmax; ymin ymax ymax ymin]';
pv2=[cos(phi); sin(phi)]';

% [p,t]=polymesh({pv1,pv2},[1,1],[1,0;0,1],[8/m,1.3]);
% [p,t] = fixmesh(p,t);

pv = [pv1; pv2];
[p,t]=polymesh({pv},1,[1,0],[8/m,1.3]);

s1 = strcat('all(p(:,2)<',int2str(eval('ymin')),'+1e-6)');
s2 = strcat('all(p(:,1)>',int2str(eval('xmax')),'-1e-6)');
s3 = strcat('all(p(:,2)>',int2str(eval('ymax')),'-1e-6)');
s4 = strcat('all(p(:,1)<',int2str(eval('xmin')),'+1e-6)');
bndexpr={s1,s2,s3,s4,'all(sum(p.^2,2)<2^2)'};   

fd=@(p) sqrt(sum(p.^2,2))-1;
mesh = mkmesh(p,t,porder,bndexpr,0,0,5,fd);


% R = 1;
% W = 5;
% L = 10;
% H = 5;
% n1= 40;
% h2= 0.5;
% growth=1.3;
% 
% phi=linspace(pi,0,n1+1)';
% phi=logdec(phi,1);
% pv1=[R*cos(phi),R*sin(phi); L,0; L,L; -L,L; -L 0];
% 
% [mesh.p,mesh.t]=polymesh({pv1},1,[1,0],[h2,growth],@hh,[]);
% [mesh.p,mesh.t] = fixmesh(mesh.p,mesh.t);
% [mesh.f,mesh.t2f] = mkt2f(mesh.t);
% 
% bndexpr={sprintf('all(sqrt(sum(p.^2,2))<%g+1e-4)',R), ...         
%          sprintf('all(p(:,2)<%g+1e-4)',0), ...
%          sprintf('all(p(:,1)>%g-1e-4)',L), ...
%          sprintf('all(p(:,2)>%g-1e-4)',L), ...          
%          sprintf('all(p(:,1)<%g+1e-4)',-L)};
%      
% % bndexpr = {'all(p(:,2)<0+1e-6)','all(p(:,1)>5-1e-6)', ...
% %            'all(p(:,2)>5-1e-6)','all(p(:,1)<-5+1e-6)','true'};     
%      
% mesh.f = setbndnbrs(mesh.p,mesh.f,bndexpr);
% mesh.fcurved = (mesh.f(:,4)==-1);
% ic = mesh.fcurved;
% mesh.tcurved = false(size(mesh.t,1),1);
% mesh.tcurved(mesh.f(ic,3)) = true;
% 
% fd=@(p) sqrt(sum(p.^2,2))-R;
% mesh.porder = porder;
% %[mesh.plocal,mesh.tlocal] = uniformlocalpnts(mesh.porder);
% %mesh.dgnodes = createnodes(mesh,fd);
% [mesh.plocal,mesh.tlocal,mesh.plocfc,mesh.tlocfc,permnode,permedge,permface] = masternodes(mesh.porder,2,0,0);
% mesh.permnode = permnode;
% mesh.permedge = permedge;
% mesh.permface = permface;
% mesh.perm = permedge;
% 
% mesh.dgnodes = createnodes(mesh,fd);
% mesh.ne = size(mesh.t,1);
% mesh.np = size(mesh.p,1);
% mesh.nf = size(mesh.f,1);
% mesh.nd = 2;
% mesh.elemtype = 0;
% mesh.nodetype = 0;
% 
% 
% function h=hh(p,hparabola)
% 
% h=inf+0*p(:,1);
% if ~isempty(hparabola)
%   p(:,1)=p(:,1)-hparabola(1);
%   p=p(:,[2,1]);
%   h1=hparabola(3)+hparabola(4)*abs(dparabola(p,hparabola(2)));
%   h=min(h,h1);
% end
