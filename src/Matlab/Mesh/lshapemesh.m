function [p,t] = lshapemesh(m,elemtype)

a=2;
n = 2*m;
[p1,t1]=squaremesh(n,m,1,elemtype);
p1(1,:)=(a)*p1(1,:);
p1(2,:)=(a/2)+(a/2)*p1(2,:);

[p2,t2]=squaremesh(m,m,1,elemtype);

[p,t] = connectmesh(p1,t1,p2,t2);
p = p - 1;
