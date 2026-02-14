function m = bnd_ellipsoid(p, param, varargin)
%BND_ELLIPSOID Predicate for an ellipsoid.
%   m = bnd_ellipsoid(p, param) returns logical row vector m for points on:
%       sum_i ((q_i / a_i)^2) = 1  (within tolerance),
%   where q = p - center for axis-aligned ellipsoids.
%
%   m = bnd_ellipsoid(p, param, mode) supports:
%       mode = "boundary" (default), "inside", "outside",
%              "strict-inside", "strict-outside"
%
%   Supported parameterizations:
%   1) Numeric:
%      [center(1:dim), axes(1:dim)] or [center(1:dim), axes(1:dim), tol]
%   2) Struct:
%      param.center   : dim-by-1 center
%      param.axes     : dim-by-1 semi-axis lengths (positive)
%      param.tol      : optional tolerance (default 1e-8)
%      param.mode     : optional mode string
%      param.rotation : optional dim-by-dim rotation matrix.
%                       If provided, q = rotation.' * (p - center).
%
%   p is dim-by-np.

dim = size(p,1);
tol = 1e-8;
mode = "boundary";
rotation = [];

if isnumeric(param)
    if numel(param) ~= 2*dim && numel(param) ~= 2*dim + 1
        error("bnd_ellipsoid: numeric param must be [center(1:dim), axes(1:dim), (tol)].");
    end
    center = param(1:dim).';
    axesv = param(dim+1:2*dim).';
    if numel(param) == 2*dim + 1
        tol = param(end);
    end
elseif isstruct(param)
    if ~isfield(param, "center") || ~isfield(param, "axes")
        error("bnd_ellipsoid: struct param requires center and axes.");
    end
    center = param.center(:);
    axesv = param.axes(:);
    if isfield(param, "tol")
        tol = param.tol;
    end
    if isfield(param, "mode")
        mode = param.mode;
    end
    if isfield(param, "rotation")
        rotation = param.rotation;
    end
else
    error("bnd_ellipsoid: param must be numeric or struct.");
end

if nargin >= 3
    mode = varargin{1};
end

if numel(center) ~= dim || numel(axesv) ~= dim
    error("bnd_ellipsoid: center/axes lengths must equal size(p,1).");
end
if any(axesv <= 0)
    error("bnd_ellipsoid: all axes lengths must be positive.");
end

q = p - center;
if ~isempty(rotation)
    if ~ismatrix(rotation) || any(size(rotation) ~= [dim dim])
        error("bnd_ellipsoid: rotation must be dim-by-dim.");
    end
    q = rotation.' * q;
end

phi = sum((q ./ axesv).^2, 1);
m = applyMode(phi, tol, mode);
end

function m = applyMode(phi, tol, mode)
if isnumeric(mode)
    s = sign(mode);
    if s == 0
        m = abs(phi - 1) <= tol;
    elseif s > 0
        m = phi <= 1 + tol;
    else
        m = phi >= 1 - tol;
    end
    return;
end

if ~(ischar(mode) || isstring(mode))
    error("bnd_ellipsoid: mode must be numeric, char, or string.");
end

ms = lower(strtrim(char(mode)));
if strcmp(ms, "boundary") || strcmp(ms, "on") || strcmp(ms, "surface") || strcmp(ms, "=")
    m = abs(phi - 1) <= tol;
elseif strcmp(ms, "inside") || strcmp(ms, "in") || strcmp(ms, "interior") || strcmp(ms, "<") || strcmp(ms, "<=")
    m = phi <= 1 + tol;
elseif strcmp(ms, "outside") || strcmp(ms, "out") || strcmp(ms, "exterior") || strcmp(ms, ">") || strcmp(ms, ">=")
    m = phi >= 1 - tol;
elseif strcmp(ms, "strict-inside") || strcmp(ms, "inside-strict")
    m = phi < 1 - tol;
elseif strcmp(ms, "strict-outside") || strcmp(ms, "outside-strict")
    m = phi > 1 + tol;
else
    error("bnd_ellipsoid: unsupported mode '%s'.", ms);
end
end
