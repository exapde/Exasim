function p = mkline(p1, p2, n)
%MKLINE Generate n equally spaced points between p1 and p2
%
%   p = mkline(p1, p2, n)
%
%   Inputs:
%       p1 : 1×d vector (starting point)
%       p2 : 1×d vector (ending point)
%       n  : number of points (>=2)
%
%   Output:
%       p  : n×d matrix of points along the line

    if nargin ~= 3
        error('mkline requires three inputs: p1, p2, n.');
    end

    if n < 2
        error('n must be at least 2.');
    end

    p1 = p1(:)';   % ensure row vector
    p2 = p2(:)';

    if length(p1) ~= length(p2)
        error('p1 and p2 must have the same dimension.');
    end

    t = linspace(0, 1, n)';   % parameter from 0 to 1
    p = p1 + t .* (p2 - p1);

end