function mesh = mkmesh_bullet(m,n,porder)

% mesh = mkmesh_bullet([16 6],[6 16 16],3);
% mesh = mkmesh_bullet([31 11],[11 31 31],3);

elemtype = 1;
nodetype = 1;

a    = 8;
b    = 12;
c    = 1;
d    = 1;

cw   = 1.2;
cr   = 1.05;

[p1,t1] = squaremesh(m(1),n(1),0,elemtype);
p1(:,1) = loginc(p1(:,1),cw);
p1(:,1) = c+b*p1(:,1);
p1(:,2) = (d/2)*logdec(p1(:,2),cw);

[p2,t2] = squaremesh(m(1),n(2),0,elemtype);
p2(:,1) = loginc(p2(:,1),cw);
p2(:,1) = c+b*p2(:,1);
p2(:,2) = exp(log(2*(a+d/2))*p2(:,2).^cr)/2;

[p3,t3] = squaremesh(m(2),n(2),0,elemtype);
p3(:,1) = c*logdec(loginc(p3(:,1),cw),cw);
p3(:,2) = exp(log(2*(a+d/2))*p3(:,2).^cr)/2;

[p4,t4] = squaremesh(n(2),n(3),0,elemtype);
bndexpr = {'all(p(:,2)<1e-3)','all(p(:,1)>max(p0(:,1))-1e-3)', ...
           'all(p(:,2)>max(p0(:,2))-1e-3)','all(p(:,1)<1e-3)'};     
mesh4 = mkmesh(p4,t4,porder,bndexpr,elemtype,nodetype);
x = log(2*(a+d/2))*p4(:,1).^cr;
y = pi/2+pi*p4(:,2);
p4 = [exp(x).*cos(y) exp(x).*sin(y)]/2;

ind = [1 4 3 2];
%ind = [1 2 3];
p5 = p3; p5(:,2) = -p5(:,2);
t5 = t3(:,ind);

p6 = p2; p6(:,2) = -p6(:,2);
t6 = t2(:,ind);

p7 = p1; p7(:,2) = -p7(:,2);
t7 = t1(:,ind);

[p,t] = connectmesh(p1,t1,p2,t2);
[p,t] = connectmesh(p,t,p3,t3);
[p,t] = connectmesh(p,t,p4,t4);
[p,t] = connectmesh(p,t,p5,t5);
[p,t] = connectmesh(p,t,p6,t6);
[p,t] = connectmesh(p,t,p7,t7);

bndexpr = {'all(sqrt(sum(p.^2,2))<3)','all(sqrt(sum(p.^2,2))>3)'};     
mesh = mkmesh(p,t,porder,bndexpr,elemtype,nodetype);

x = mesh4.dgnodes(:,1,:);
y = mesh4.dgnodes(:,2,:);
x = log(2*(a+d/2))*x.^cr;
y = pi/2+pi*y;
it = size(t1,1)+size(t2,1)+size(t3,1);
mesh.dgnodes(:,1,it+1:it+size(t4,1)) = exp(x).*cos(y)/2;
mesh.dgnodes(:,2,it+1:it+size(t4,1)) = exp(x).*sin(y)/2;

   