function m = bnd_or(varargin)
%BND_OR Logical OR for boundary predicates.
%   m = bnd_or(a, b, c, ...) where each input is logical-compatible.

if nargin == 0
    m = false;
    return;
end

m = false(size(varargin{1}));
for k = 1:nargin
    m = m | logical(varargin{k});
end
end
