function [sr,sr_udg] = source_ns(p,udg,param,time)

[ng,nc] = size(udg);
nch = 4;

sr = zeros(ng,nch);
sr_udg = zeros(ng,nch,nc); 
