function m = bnd_halfspace(p, param)
%BND_HALFSPACE Boundary predicate for a half-space.
%   m = bnd_halfspace(p, param) returns logical row vector m where points
%   satisfy a linear inequality.
%
%   Supported parameterizations:
%   1) Axis-aligned (numeric): [axis, value, side] or [axis, value, side, tol]
%      side = +1 for >=, side = -1 for <=
%   2) Axis-aligned (struct):  param.axis, param.value, param.side, optional tol
%   3) General (struct):       param.normal with either param.point or param.offset,
%      plus param.side and optional param.tol
%
%   p is dim-by-np.

tol = 0;

if isnumeric(param)
    if numel(param) < 3
        error("bnd_halfspace: numeric param must be [axis, value, side, (tol)].");
    end
    axis = round(param(1));
    value = param(2);
    side = parseSide(param(3));
    if numel(param) >= 4
        tol = param(4);
    end
    checkAxis(axis, size(p,1));
    d = p(axis,:) - value;
    m = applySide(d, side, tol);
    return;
end

if ~isstruct(param)
    error("bnd_halfspace: param must be numeric or struct.");
end

if isfield(param, "tol")
    tol = param.tol;
end

if isfield(param, "axis") && isfield(param, "value")
    axis = round(param.axis);
    value = param.value;
    checkAxis(axis, size(p,1));
    side = parseSide(param.side);
    d = p(axis,:) - value;
    m = applySide(d, side, tol);
    return;
end

if isfield(param, "normal")
    normal = param.normal(:);
    if numel(normal) ~= size(p,1)
        error("bnd_halfspace: normal length must equal size(p,1).");
    end

    side = parseSide(param.side);
    if isfield(param, "point")
        point = param.point(:);
        if numel(point) ~= size(p,1)
            error("bnd_halfspace: point length must equal size(p,1).");
        end
        d = normal.' * (p - point);
    elseif isfield(param, "offset")
        d = normal.' * p + param.offset;
    else
        error("bnd_halfspace: struct with normal requires point or offset.");
    end

    m = applySide(d, side, tol);
    return;
end

error("bnd_halfspace: unsupported param format.");
end

function m = applySide(d, side, tol)
if side >= 0
    m = d >= -tol;
else
    m = d <= tol;
end
end

function side = parseSide(sideValue)
if isnumeric(sideValue)
    side = sign(sideValue);
    if side == 0
        error("bnd_halfspace: side must be non-zero (+1 or -1).");
    end
    return;
end

if isstring(sideValue) || ischar(sideValue)
    s = char(sideValue);
    if strcmp(s, ">=") || strcmp(s, ">")
        side = 1;
        return;
    end
    if strcmp(s, "<=") || strcmp(s, "<")
        side = -1;
        return;
    end
end

error("bnd_halfspace: side must be +1/-1 or one of '>=', '>', '<=', '<'.");
end

function checkAxis(axis, dim)
if axis < 1 || axis > dim
    error("bnd_halfspace: axis must be between 1 and size(p,1).");
end
end
