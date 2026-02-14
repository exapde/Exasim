function m = bnd_le(p, param)
%BND_LE Coordinate less-than-or-equal predicate.
%   m = bnd_le(p, param), param = [axis, value] or [axis, value, tol]
%   Also supports struct with fields axis, value, optional tol.

[axis, value, tol] = parseAxisParam(p, param, "bnd_le");
m = p(axis,:) <= value + tol;
end

function [axis, value, tol] = parseAxisParam(p, param, fname)
tol = 0;
if isnumeric(param)
    if numel(param) < 2
        error("%s: numeric param must be [axis, value] or [axis, value, tol].", fname);
    end
    axis = round(param(1));
    value = param(2);
    if numel(param) >= 3
        tol = param(3);
    end
elseif isstruct(param)
    if ~isfield(param, "axis") || ~isfield(param, "value")
        error("%s: struct param requires axis and value.", fname);
    end
    axis = round(param.axis);
    value = param.value;
    if isfield(param, "tol")
        tol = param.tol;
    end
else
    error("%s: param must be numeric or struct.", fname);
end

if axis < 1 || axis > size(p,1)
    error("%s: axis must be between 1 and size(p,1).", fname);
end
end
