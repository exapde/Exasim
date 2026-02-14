function m = bnd_plane(p, param, varargin)
%BND_PLANE Predicate for a plane with side selection.
%   m = bnd_plane(p, param) returns points on the plane (default mode).
%
%   m = bnd_plane(p, param, mode) supports:
%       mode = "boundary" (default), "plus", "minus",
%              "strict-plus", "strict-minus"
%
%   Supported parameterizations:
%   1) Axis-aligned (numeric): [axis, value] or [axis, value, tol]
%   2) Axis-aligned (struct):  param.axis, param.value, optional param.tol
%   3) General (struct):       param.normal with either
%      - param.point  (point on plane), or
%      - param.offset (normal.' * x + offset = 0)
%      Optional: param.tol, param.mode
%
%   p is dim-by-np.

tol = 1e-8;
mode = "boundary";

if isnumeric(param)
    if numel(param) < 2
        error("bnd_plane: numeric param must be [axis, value] or [axis, value, tol].");
    end
    axis = round(param(1));
    value = param(2);
    if numel(param) >= 3
        tol = param(3);
    end
    checkAxis(axis, size(p,1));
    d = p(axis,:) - value;
    if nargin >= 3
        mode = varargin{1};
    end
    m = applyMode(d, tol, mode);
    return;
end

if ~isstruct(param)
    error("bnd_plane: param must be numeric or struct.");
end

if isfield(param, "tol")
    tol = param.tol;
end
if isfield(param, "mode")
    mode = param.mode;
end
if nargin >= 3
    mode = varargin{1};
end

if isfield(param, "axis") && isfield(param, "value")
    axis = round(param.axis);
    value = param.value;
    checkAxis(axis, size(p,1));
    d = p(axis,:) - value;
    m = applyMode(d, tol, mode);
    return;
end

if isfield(param, "normal")
    normal = param.normal(:);
    if numel(normal) ~= size(p,1)
        error("bnd_plane: normal length must equal size(p,1).");
    end

    if isfield(param, "point")
        point = param.point(:);
        if numel(point) ~= size(p,1)
            error("bnd_plane: point length must equal size(p,1).");
        end
        d = normal.' * (p - point);
        m = applyMode(d, tol, mode);
        return;
    end

    if isfield(param, "offset")
        d = normal.' * p + param.offset;
        m = applyMode(d, tol, mode);
        return;
    end

    error("bnd_plane: struct with normal requires either point or offset.");
end

error("bnd_plane: unsupported param format.");
end

function checkAxis(axis, dim)
if axis < 1 || axis > dim
    error("bnd_plane: axis must be between 1 and size(p,1).");
end
end

function m = applyMode(d, tol, mode)
if isnumeric(mode)
    s = sign(mode);
    if s == 0
        m = abs(d) <= tol;
    elseif s > 0
        m = d >= -tol;
    else
        m = d <= tol;
    end
    return;
end

if ~(ischar(mode) || isstring(mode))
    error("bnd_plane: mode must be numeric, char, or string.");
end

ms = lower(strtrim(char(mode)));
if strcmp(ms, "boundary") || strcmp(ms, "on") || strcmp(ms, "surface") || strcmp(ms, "=")
    m = abs(d) <= tol;
elseif strcmp(ms, "plus") || strcmp(ms, "+") || strcmp(ms, "positive") || strcmp(ms, ">") || strcmp(ms, ">=")
    m = d >= -tol;
elseif strcmp(ms, "minus") || strcmp(ms, "-") || strcmp(ms, "negative") || strcmp(ms, "<") || strcmp(ms, "<=")
    m = d <= tol;
elseif strcmp(ms, "strict-plus") || strcmp(ms, "plus-strict")
    m = d > tol;
elseif strcmp(ms, "strict-minus") || strcmp(ms, "minus-strict")
    m = d < -tol;
else
    error("bnd_plane: unsupported mode '%s'.", ms);
end
end
