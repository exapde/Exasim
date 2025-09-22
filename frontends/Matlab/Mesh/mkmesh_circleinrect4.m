function mesh=mkmesh_circleinrect4(porder,mul1,mul2)

if nargin<2, mul1=1; end
if nargin<3, mul2=1; end

R = 0.5;
h = 0.35;

xmin = -10;
xmax =  10;
ymin = -3;
ymax =  3;

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
[p,t]=polymesh({pv1,pv2,pv3},[1,1,1],[1,0;0,1;0,1],[h,1.5]);
[p,t] = fixmesh(p,t);

s1 = strcat('all(p(:,2)<',int2str(eval('ymin')),'+1e-6)');
s2 = strcat('all(p(:,1)>',int2str(eval('xmax')),'-1e-6)');
s3 = strcat('all(p(:,2)>',int2str(eval('ymax')),'-1e-6)');
s4 = strcat('all(p(:,1)<',int2str(eval('xmin')),'+1e-6)');
bndexpr={s1,s2,s3,s4,'all(sum(p.^2,2)<4^2)'};   

%fd=@(p) sqrt(sum(p.^2,2))-R;
fd=@(p) min(sqrt(p(:,1).^2 + (p(:,2)-1.5*R).^2)-R,sqrt(p(:,1).^2 + (p(:,2)+1.5*R).^2)-R);
mesh = mkmesh(p,t,porder,bndexpr,0,1,5,fd);
