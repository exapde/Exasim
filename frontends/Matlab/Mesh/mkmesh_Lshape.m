function mesh = mkmesh_Lshape(porder,m,elemtype,nodetype)

% 
% a=2;
% n = 2*m-1;
% [p1,t1]=squaremesh(m,n);
% p1(:,1)=(a/2)*p1(:,1);
% p1(:,2)=a*p1(:,2);
% 
% [p2,t2]=squaremesh(m,m);
% p2(:,1)=a/2+(a/2)*p2(:,1);
% p2(:,2)=(a/2)*p2(:,2);
% 
% [p,t] = connectmesh(p1,t1,p2,t2);
% p(:,1) = p(:,1)-1;
% p(:,2) = p(:,2)-1;

a=2;
n = 2*m-1;
[p1,t1]=squaremesh(n,m,1,elemtype);
p1(:,1)=(a)*p1(:,1);
p1(:,2)=(a/2)+(a/2)*p1(:,2);

[p2,t2]=squaremesh(m,m,1,elemtype);

[p,t] = connectmesh(p1,t1,p2,t2);
p(:,1) = p(:,1)-1;
p(:,2) = p(:,2)-1;

bndexpr = {'true'};     
mesh = mkmesh(p,t,porder,bndexpr,elemtype,nodetype);


% a=2;
% n = 2*m-1;
% [p1,t1]=squaremesh(n,m);
% p1(:,1)=(a)*p1(:,1);
% p1(:,2)=(a/2)+(a/2)*p1(:,2);
% 
% [p2,t2]=squaremesh(m,m);
% 
% [p,t] = connectmesh(p1,t1,p2,t2);
% p(:,1) = p(:,1)-1;
% p(:,2) = p(:,2)-1;