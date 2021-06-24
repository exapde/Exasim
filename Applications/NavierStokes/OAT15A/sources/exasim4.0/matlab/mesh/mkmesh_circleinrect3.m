function mesh=mkmesh_circleinrect3(porder,mul1,mul2)

if nargin<2, mul1=1; end
if nargin<3, mul2=1; end

hs=[.15,.15,.15,0.25,0.5,2];
R = 0.5;

xmin = -20/2;
xmax = 60/2;
ymin = -24/2;
ymax = 24/2;

m=4*mul1;
n=12*mul1;
o=4*mul2;

s1=(1:m)/m; z1=0*s1;
s2=(1:n)/n; z2=0*s2;
phi=2*pi*(1:4*o)/(4*o);

pv1=[xmin+(xmax-xmin)*s2,xmax+z1,xmax-(xmax-xmin)*s2,xmin+z1;
     ymin+z2,ymin+(ymax-ymin)*s1,ymax+z2,ymax-(ymax-ymin)*s1]';
pva=R*[cos(phi); sin(phi)]';
pv2 = pva; 
pv2(:,2) = pva(:,2)+1.5*R;
pv3 = pva; 
pv3(:,2) = pva(:,2)-1.5*R;

%[p,t]=polymesh({pv1,pv2},[1,1],[1,0;0,1],[8/m,1.3]);
[p,t]=polymesh({pv1,pv2,pv3},[1,1,1],[1,0;0,1;0,1],[hs(end),1.5],@href,hs);
[p,t] = fixmesh(p,t);

s1 = strcat('all(p(:,2)<',int2str(eval('ymin')),'+1e-6)');
s2 = strcat('all(p(:,1)>',int2str(eval('xmax')),'-1e-6)');
s3 = strcat('all(p(:,2)>',int2str(eval('ymax')),'-1e-6)');
s4 = strcat('all(p(:,1)<',int2str(eval('xmin')),'+1e-6)');
bndexpr={s1,s2,s3,s4,'all(sum(p.^2,2)<4^2)'};   

%fd=@(p) sqrt(sum(p.^2,2))-R;
fd=@(p) min(sqrt(p(:,1).^2 + (p(:,2)-1.5*R).^2)-R,sqrt(p(:,1).^2 + (p(:,2)+1.5*R).^2)-R);
mesh = mkmesh(p,t,porder,bndexpr,0,1,5,fd);

function h=href(p,hs)
% c = 5;
% pv1=[20,c;-0.5,2.5;-0.5,-2;20,-c;20,c];
% pv2=[5,c;-0.5,2.5;-0.5,-2;5,-c;5,c];
pv0=[30,9;20,7;20,-7;30,-9;30,9];
pv1=[20,7;10,4.5;10,-4.5;20,-7;20,7];
pv2=[10,4.5;3,3;3,-3;10,-4.5;10,4.5];
pv3=[3,3;-0.5,1.5;-0.5,-1.5;3,-3;3,3];
inside0=inpolygon(p(:,1),p(:,2),pv0(:,1),pv0(:,2));
inside1=inpolygon(p(:,1),p(:,2),pv1(:,1),pv1(:,2));
inside2=inpolygon(p(:,1),p(:,2),pv2(:,1),pv2(:,2));
inside3=inpolygon(p(:,1),p(:,2),pv3(:,1),pv3(:,2));
h=inf*ones(size(p,1),1);
h(inside0)=1.5*hs(5);
h(inside1)=hs(5);
h(inside2)=0.5*(hs(4)+hs(5));
h(inside3)=0.8*hs(4);

