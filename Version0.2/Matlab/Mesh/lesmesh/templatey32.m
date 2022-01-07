function [p,t] = templatey32(x1, y1, z1, x2, y2, z2, y3)

[p3,t3] = templatey3(x1, y1, z1, x2, y2, z2);
[p2,t2] = templatey2(x1, y2, z1, x2, y3, z2);

p = [p3; p2];
t = [t3; t2+size(p3,1)];

[p,t]=fixmesh2(p,t);

% ptplot(p,t);
% xlabel('x');
% ylabel('y');
% zlabel('z');
% axis on;

