function [sr,sr_udg] = source(p,udg,param,time)

[ng,nc] = size(udg);
nch = 1;

% sr = zeros(ng,nch);
    
xg = p(:,1);
yg = p(:,2);
zg = p(:,3);
sr = 3*pi^2*sin(pi*xg).*sin(pi*yg).*sin(pi*zg);
sr_udg = zeros(ng,nch,nc); 
