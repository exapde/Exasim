function [sr,sr_udg] = source(p,udg,param,time)

[ng,nc] = size(udg);
nch = 1;

% sr = zeros(ng,nch);
    
xg = p(:,1);
yg = p(:,2);
sr = xg.*sin(yg);
sr_udg = zeros(ng,nch,nc); 
