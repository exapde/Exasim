function [sr,sr_udg] = source(p,udg,param,time)

[ng,nc] = size(udg);
nch = 1;

% sr = zeros(ng,nch);
    
xg = p(:,1);
yg = p(:,2);
sr = 2*pi^2*sin(pi*xg).*sin(pi*yg);
sr_udg = zeros(ng,nch,nc); 

