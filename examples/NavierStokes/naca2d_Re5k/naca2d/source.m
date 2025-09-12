function [sr,sr_udg] = source(p,udg,param,time)

[ng,nc] = size(udg);

if nc==12
    nch = 4;
elseif nc==20
    nch = 5;
end

sr = zeros(ng,nch);

sr_udg = zeros(ng,nch,nc); 
