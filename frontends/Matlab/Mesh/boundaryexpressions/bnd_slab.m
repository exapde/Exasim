function m = bnd_slab(p, param, varargin)
%BND_SLAB Predicate for a slab between two parallel planes.
%   m = bnd_slab(p, param) returns points on the slab boundary (default).
%
%   m = bnd_slab(p, param, mode) supports:
%       mode = "boundary" (default), "inside", "outside",
%              "strict-inside", "strict-outside"
%
%   Supported parameterizations:
%   1) Numeric (axis-aligned): [axis, min, max] or [axis, min, max, tol]
%   2) Struct axis-aligned:    param.axis, param.min, param.max
%   3) Struct general:
%      - param.normal, param.point, param.tmin, param.tmax
%        where t = normal_unit.' * (p - point)
%      - or param.normal, param.dmin, param.dmax
%        where d = normal_unit.' * p
%
%   Optional struct fields: param.tol, param.mode.
%   p is dim-by-np.

dim = size(p,1);
tol = 1e-8;
mode = "boundary";

if isnumeric(param)
    if numel(param) < 3
        error("bnd_slab: numeric param must be [axis, min, max, (tol)].");
    end
    axis = round(param(1));
    checkAxis(axis, dim);
    d = p(axis,:);
    dmin = param(2);
    dmax = param(3);
    if numel(param) >= 4
        tol = param(4);
    end
else
    if ~isstruct(param)
        error("bnd_slab: param must be numeric or struct.");
    end
    if isfield(param, "tol")
        tol = param.tol;
    end
    if isfield(param, "mode")
        mode = param.mode;
    end

    if isfield(param, "axis") && isfield(param, "min") && isfield(param, "max")
        axis = round(param.axis);
        checkAxis(axis, dim);
        d = p(axis,:);
        dmin = param.min;
        dmax = param.max;
    elseif isfield(param, "normal")
        n = param.normal(:);
        if numel(n) ~= dim
            error("bnd_slab: normal length must equal size(p,1).");
        end
        nn = norm(n);
        if nn == 0
            error("bnd_slab: normal must be non-zero.");
        end
        nu = n / nn;

        if isfield(param, "point") && isfield(param, "tmin") && isfield(param, "tmax")
            q0 = param.point(:);
            if numel(q0) ~= dim
                error("bnd_slab: point length must equal size(p,1).");
            end
            d = nu.' * (p - q0);
            dmin = param.tmin;
            dmax = param.tmax;
        elseif isfield(param, "dmin") && isfield(param, "dmax")
            d = nu.' * p;
            dmin = param.dmin;
            dmax = param.dmax;
        else
            error("bnd_slab: with normal, provide (point,tmin,tmax) or (dmin,dmax).");
        end
    else
        error("bnd_slab: unsupported struct format.");
    end
end

if nargin >= 3
    mode = varargin{1};
end

if dmax <= dmin
    error("bnd_slab: require max > min.");
end

insideLoose = (d >= dmin - tol) & (d <= dmax + tol);
insideStrict = (d > dmin + tol) & (d < dmax - tol);
onBoundary = insideLoose & ((abs(d - dmin) <= tol) | (abs(d - dmax) <= tol));
outsideLoose = ~insideStrict;
outsideStrict = (d < dmin - tol) | (d > dmax + tol);

m = selectMode(onBoundary, insideLoose, insideStrict, outsideLoose, outsideStrict, mode);
end

function m = selectMode(onBoundary, insideLoose, insideStrict, outsideLoose, outsideStrict, mode)
if isnumeric(mode)
    s = sign(mode);
    if s == 0
        m = onBoundary;
    elseif s > 0
        m = insideLoose;
    else
        m = outsideLoose;
    end
    return;
end

if ~(ischar(mode) || isstring(mode))
    error("bnd_slab: mode must be numeric, char, or string.");
end

ms = lower(strtrim(char(mode)));
if strcmp(ms, "boundary") || strcmp(ms, "on") || strcmp(ms, "surface") || strcmp(ms, "=")
    m = onBoundary;
elseif strcmp(ms, "inside") || strcmp(ms, "in") || strcmp(ms, "interior") || strcmp(ms, "<") || strcmp(ms, "<=")
    m = insideLoose;
elseif strcmp(ms, "outside") || strcmp(ms, "out") || strcmp(ms, "exterior") || strcmp(ms, ">") || strcmp(ms, ">=")
    m = outsideLoose;
elseif strcmp(ms, "strict-inside") || strcmp(ms, "inside-strict")
    m = insideStrict;
elseif strcmp(ms, "strict-outside") || strcmp(ms, "outside-strict")
    m = outsideStrict;
else
    error("bnd_slab: unsupported mode '%s'.", ms);
end
end

function checkAxis(axis, dim)
if axis < 1 || axis > dim
    error("bnd_slab: axis must be between 1 and size(p,1).");
end
end
