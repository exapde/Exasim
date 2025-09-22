function mesh = mkmesh_freearc2daxis_aircraftskin(porder,n,a)

% Example: mesh = mkmesh_freearc2daxis(2,[11 51 51 11 15],[2.25 2.5 1.5]);

% geometry parameters
Lref = 1000;
L = 50/Lref; % r width  
H = 50/Lref; % z width
G = 2/Lref;  % aircraft skin depth
A = 7/Lref;  % air extension
C = 2/Lref;  % plasma column width
D = H-A-G;   % plasma column depth
B = L-C;

% mesh resolution parameters
n1 = n(1); 
n2 = n(2); 
n3 = n(3);
n4 = n(4); 
n5 = n(5); 

% scaling factors
a1 = a(1); 
a2 = a(2);
a3 = a(3);
 
% quad elements
elemtype = 1;
nodetype = 1;

% plasma column
[p1,t1] = squaremesh(n1,n2,0,elemtype);
p1(:,2) = loginc(p1(:,2),a1);
p1(:,1) = C*p1(:,1);
p1(:,2) = D*p1(:,2);

% upper air
[p2,t2] = squaremesh(n3,n2,0,elemtype);
p2(:,1) = loginc(p2(:,1),a2);
p2(:,2) = loginc(p2(:,2),a1);
p2(:,1) = C+B*p2(:,1);
p2(:,2) = D*p2(:,2);

% aircraft skin under air
[p3,t3] = squaremesh(n3,n4,0,elemtype);
p3(:,1) = loginc(p3(:,1),a2);
p3(:,1) = C+B*p3(:,1);
p3(:,2) = -G+G*p3(:,2);

% aircraft skin under plasma column
[p4,t4] = squaremesh(n1,n4,0,elemtype);
p4(:,1) = C*p4(:,1);
p4(:,2) = -G+G*p4(:,2);

% air under aircfrat skin

[p5,t5] = squaremesh(n3,n5,0,elemtype);
p5(:,1) = loginc(p5(:,1),a2);
p5(:,2) = logdec(p5(:,2),a3);
p5(:,1) = C+B*p5(:,1);
p5(:,2) = -(G+A)+A*p5(:,2);

% air under aircfrat skin
[p6,t6] = squaremesh(n1,n5,0,elemtype);
p6(:,2) = logdec(p6(:,2),a3);
p6(:,1) = C*p6(:,1);
p6(:,2) = -(G+A)+A*p6(:,2);

[p,t] = connectmesh(p1,t1,p2,t2);
[p,t] = connectmesh(p,t,p3,t3);
[p,t] = connectmesh(p,t,p4,t4);
[p,t] = connectmesh(p,t,p5,t5);
[p,t] = connectmesh(p,t,p6,t6);

bndexpr = {'all(p(:,2)<min(p0(:,2))+1e-4)',...
           'all(p(:,1)>max(p0(:,1))-1e-4)',...
           'all(p(:,2)>max(p0(:,2))-1e-4)',...
           'all(p(:,1)<min(p0(:,1))+1e-4)',...
           'true'};     
mesh = mkmesh(p,t,porder,bndexpr,elemtype,nodetype);

close all
meshplot(mesh);

return;
