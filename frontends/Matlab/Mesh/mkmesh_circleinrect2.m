function mesh=mkmesh_circleinrect2(porder,mul1,mul2)

if nargin<2, mul1=1; end
if nargin<3, mul2=1; end

hs=[.1,.1,.1,.45,.45,2];
R = 0.5;

xmin = -20/2;
xmax = 40/2;
ymin = -20/2;
ymax = 20/2;

m=4*mul1;
n=12*mul1;
o=4*mul2;

s1=(1:m)/m; z1=0*s1;
s2=(1:n)/n; z2=0*s2;
phi=2*pi*(1:4*o)/(4*o);

pv1=[xmin+(xmax-xmin)*s2,xmax+z1,xmax-(xmax-xmin)*s2,xmin+z1;
     ymin+z2,ymin+(ymax-ymin)*s1,ymax+z2,ymax-(ymax-ymin)*s1]';
pv2=R*[cos(phi); sin(phi)]';

%[p,t]=polymesh({pv1,pv2},[1,1],[1,0;0,1],[8/m,1.3]);
[p,t]=polymesh({pv1,pv2},[1,1],[1,0;0,1],[hs(end),1.5],@href,hs);
[p,t] = fixmesh(p,t);

s1 = strcat('all(p(:,2)<',int2str(eval('ymin')),'+1e-6)');
s2 = strcat('all(p(:,1)>',int2str(eval('xmax')),'-1e-6)');
s3 = strcat('all(p(:,2)>',int2str(eval('ymax')),'-1e-6)');
s4 = strcat('all(p(:,1)<',int2str(eval('xmin')),'+1e-6)');
bndexpr={s1,s2,s3,s4,'all(sum(p.^2,2)<2^2)'};   

fd=@(p) sqrt(sum(p.^2,2))-R;
mesh = mkmesh(p,t,porder,bndexpr,0,1,5,fd);

function h=href(p,hs)
pv1=[20,2.5;0,2;0,-2;20,-2.5;20,2.5];
%pv2=[10,2.5;0,2;0,-2;10,-2.5;10,2.5];
inside1=inpolygon(p(:,1),p(:,2),pv1(:,1),pv1(:,2));
%inside2=inpolygon(p(:,1),p(:,2),pv2(:,1),pv2(:,2));
h=inf*ones(size(p,1),1);
h(inside1)=hs(5);
%h(inside2)=hs(4);
