function [mesh1, mesh2] = rightleft2d(mesh1, mesh2)

a = 0.5*(mesh1.cgright + mesh2.cgleft)';
mesh1.p(:,mesh1.cright) = a;
mesh2.p(:,mesh2.cleft) = a;

a = 0.5*(mesh1.dgright + mesh2.dgleft);
x = mesh1.dgnodes(:,1,:);
y = mesh1.dgnodes(:,2,:);
x(mesh1.dright) = a(:,1);
y(mesh1.dright) = a(:,2);
mesh1.dgnodes(:,1,:) = x;
mesh1.dgnodes(:,2,:) = y;

x = mesh2.dgnodes(:,1,:);
y = mesh2.dgnodes(:,2,:);
x(mesh2.dleft) = a(:,1);
y(mesh2.dleft) = a(:,2);
mesh2.dgnodes(:,1,:) = x;
mesh2.dgnodes(:,2,:) = y;

