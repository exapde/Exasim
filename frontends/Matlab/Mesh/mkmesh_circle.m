function mesh = mkmesh_circle(porder,elemtype,nodetype,R0,H0,n1,n2)

% porder = 1;
% elemtype = 1;
% nodetype = 1;

% R0 = 1; 
% H0 = R0/3;
%n1 = 11; % number of points in the inner square
%n2 = 9;  % number of points in the radius direction 

%[p0,t0]=squaremesh(n1,n2,1,elemtype);
mesh = mkmesh_square(n1,n2,porder,1,1,1,elemtype,nodetype);
p0 = mesh.p;
t0 = mesh.t;
xi = p0(:,1);
et = p0(:,2);
theta = 3*pi/4 - pi*xi/2;
r = et.*(R0 - H0./sin(theta)) + H0./sin(theta);
x = r.*cos(theta);
y = r.*sin(theta);

p1 = [x y];
p2 = [cos(pi/2)*x-sin(pi/2)*y sin(pi/2)*x+cos(pi/2)*y];
p3 = [cos(pi)*x-sin(pi)*y sin(pi)*x+cos(pi)*y];
p4 = [cos(3*pi/2)*x-sin(3*pi/2)*y sin(3*pi/2)*x+cos(3*pi/2)*y];
t1 = t0; t2 = t0; t3 = t0; t4 = t0;

Xi = mesh.dgnodes(:,1,:);
Et = mesh.dgnodes(:,2,:);
theta = 3*pi/4 - pi*Xi/2;
r = Et.*(R0 - H0./sin(theta)) + H0./sin(theta);
X = r.*cos(theta);
Y = r.*sin(theta);
dg1 = cat(2,X,Y);
dg2 = cat(2,cos(pi/2)*X-sin(pi/2)*Y,sin(pi/2)*X+cos(pi/2)*Y);
dg3 = cat(2,cos(pi)*X-sin(pi)*Y,sin(pi)*X+cos(pi)*Y);
dg4 = cat(2,cos(3*pi/2)*X-sin(3*pi/2)*Y,sin(3*pi/2)*X+cos(3*pi/2)*Y);

ind = find(abs(y-H0)<1e-10);
[x5,y5]=ndgrid(sort(x(ind)),sort(x(ind)));
p5 = [x5(:) y5(:)];
t5 = [1 2 n1+2 n1+1];
t5 = kron(t5,ones(n1-1,1))+kron(ones(size(t5)),(0:n1-2)'*n1);
t5 = kron(t5,ones(n1-1,1))+kron(ones(size(t5)),(0:n1-2)');
dg5 = createdgnodes(p5,t5,mesh.plocal);

[p,t] = connectmesh(p1,t1,p2,t2);
[p,t] = connectmesh(p,t,p3,t3);
[p,t] = connectmesh(p,t,p4,t4);
[p,t] = connectmesh(p,t,p5,t5);

figure(1); clf; 
hold on;
simpplot(p1,t1); 
simpplot(p2,t2); 
simpplot(p3,t3); 
simpplot(p4,t4); 
simpplot(p5,t5); 
axis on;
axis tight;

figure(2); clf; 
simpplot(p,t); 
axis on;
axis tight;

bndexpr = {'true'};     
mesh = mkmesh(p,t,porder,bndexpr,elemtype,nodetype);
mesh.dgnodes = cat(3,dg1,dg2);
mesh.dgnodes = cat(3,mesh.dgnodes,dg3);
mesh.dgnodes = cat(3,mesh.dgnodes,dg4);
mesh.dgnodes = cat(3,mesh.dgnodes,dg5);



