function mesh = mkmesh_freearc2daxis(porder,n,a)

% Example: mesh = mkmesh_freearc2daxis(2,[11 51 51 11 15],[2.25 2.5 1.5]);

% geometry parameters
Lref = 1000;
L = 50/Lref; % r width  
H = 50/Lref; % z width
C = 2/Lref;  % plasma column width

% mesh resolution parameters
n1 = n(1); 
n2 = n(2); 
n3 = n(3);

% scaling factors
a1 = a(1); 
a2 = a(2);
 
% quad elements
elemtype = 1;
nodetype = 1;

% plasma column
[p1,t1] = squaremesh(n1,n2,0,elemtype);
p1(:,2) = loginc(p1(:,2),a1);
p1(:,1) = C*p1(:,1);
p1(:,2) = H*p1(:,2);

% upper air
[p2,t2] = squaremesh(n3,n2,0,elemtype);
p2(:,1) = loginc(p2(:,1),a2);
p2(:,2) = loginc(p2(:,2),a1);
p2(:,1) = C+(L-C)*p2(:,1);
p2(:,2) = H*p2(:,2);

[p,t] = connectmesh(p1,t1,p2,t2);

bndexpr = {'all(p(:,2)<min(p0(:,2))+1e-6)',...
           'all(p(:,1)>max(p0(:,1))-1e-6)',...
           'all(p(:,2)>max(p0(:,2))-1e-6 & p(:,1)>0.002-1e-6)',...
           'all(p(:,1)<min(p0(:,1))+1e-6)',...
           'true'};     
mesh = mkmesh(p,t,porder,bndexpr,elemtype,nodetype);

% close all
% meshplot(mesh);

return;
