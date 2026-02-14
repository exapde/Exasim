function m = bnd_and(varargin)
%BND_AND Logical AND for boundary predicates.
%   m = bnd_and(a, b, c, ...) where each input is logical-compatible.

if nargin == 0
    m = true;
    return;
end

m = true(size(varargin{1}));
for k = 1:nargin
    m = m & logical(varargin{k});
end
end
