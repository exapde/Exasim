function [p,t] = templatey42(x1, y1, z1, x2, y2, z2, y3)

[p4,t4] = templatey4(x1, y1, z1, x2, y2, z2);
[p2,t2] = templatey2(x1, y2, z1, x2, y3, z2);

p = [p4; p2];
t = [t4; t2+size(p4,1)];

[p,t]=fixmesh2(p,t);

% ptplot(p,t);
% xlabel('x');
% ylabel('y');
% zlabel('z');
% axis on;

