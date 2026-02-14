function [p,t] = circlemesh(n)

[p1,t1] = quartercirclemesh(n);
p2 = p1; p2(2,:) = -p2(2,:);
p3 = p1; p3(1,:) = -p3(1,:);
p4 = p2; p4(1,:) = -p4(1,:);

[p,t] = connectmesh(p1',t1',p2',t1');
[p,t] = connectmesh(p,t,p3',t1');
[p,t] = connectmesh(p,t,p4',t1');
p = p';
t = t';

% figure(1); clf; 
% simpplot(p1,t1);
% hold on;
% simpplot(p2,t1);
% simpplot(p3,t1);
% simpplot(p4,t1);
% 
% figure(2); clf; 
% simpplot(p,t);
