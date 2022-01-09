function [p,t]=linemesh(m)
%LINEMESH 1-D Regular grid for the unit interval
%
%      P:         Node positions 
%      T:         Element indices
%      m:         Number of elements
%   Example:
%      [p,t]=linemesh(10);
%

if nargin<1, m=10; end

m = m+1;
p = linspace(0,1,m);
t = [1:1:(m-1); 2:1:m];


