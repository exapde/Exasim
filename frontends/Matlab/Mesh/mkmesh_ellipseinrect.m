function mesh=mkmesh_ellipseinrect(porder,mul1,mul2)

if nargin<2, mul1=1; end
if nargin<3, mul2=1; end

xmin = -9;
xmax = 5;
ymin = -5;
ymax = 7;
a = 0.2;
b = 0.1;

m=4*mul1/2;
n=4*mul1;
o=4*mul2;

s1=(1:m)/m; z1=0*s1;
s2=(1:n)/n; z2=0*s2;
phi=2*pi*(1:4*o)/(4*o);

pv1=[xmin+(xmax-xmin)*s2,xmax+z1,xmax-(xmax-xmin)*s2,xmin+z1;
     ymin+z2,ymin+(ymax-ymin)*s1,ymax+z2,ymax-(ymax-ymin)*s1]';
pv2=[a*cos(phi); b*sin(phi)]';

[p,t]=polymesh({pv1,pv2},[1,1],[1,0;0,1],[1.1*(xmax-xmin)/n,1.9]);
[p,t] = fixmesh(p,t);
figure(1); clf; simpplot(p,t); 
hold on; plot(pv1(:,1),pv1(:,2),'*'); 
hold on; plot(pv2(:,1),pv2(:,2),'*'); 
axis on;

s1 = strcat('all(p(:,2)<',num2str(eval('ymin')),'+1e-6)');
s2 = strcat('all(p(:,1)>',num2str(eval('xmax')),'-1e-6)');
s3 = strcat('all(p(:,2)>',num2str(eval('ymax')),'-1e-6)');
s4 = strcat('all(p(:,1)<',num2str(eval('xmin')),'+1e-6)');
bndexpr={s1,s2,s3,s4,'all(sum(p.^2,2)<1^2)'};

fd=@(p) sqrt((p(:,1)/a).^2+(p(:,2)/b).^2)-1;
mesh = mkmesh(p,t,porder,bndexpr,0,1,5,fd);

