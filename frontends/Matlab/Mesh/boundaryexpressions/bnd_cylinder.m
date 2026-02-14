function m = bnd_cylinder(p, param, varargin)
%BND_CYLINDER Predicate for an infinite or finite cylinder.
%   m = bnd_cylinder(p, param) returns logical row vector m for points on:
%       distance(point, axis line) = radius (default mode, within tolerance).
%
%   m = bnd_cylinder(p, param, mode) supports:
%       mode = "boundary" (default), "inside", "outside",
%              "strict-inside", "strict-outside"
%
%   Required struct fields:
%     param.center : point on axis (dim-by-1)
%     param.axis   : axis direction (dim-by-1), non-zero
%     param.radius : cylinder radius
%   Optional:
%     param.tol    : default 1e-8
%     param.mode   : selection mode (see above)
%     param.length : finite cylinder length centered at param.center
%     param.height : alias for param.length
%     param.tmin   : lower axial limit (relative to center along axis)
%     param.tmax   : upper axial limit (relative to center along axis)
%                    If tmin/tmax are omitted, cylinder is infinite.
%
%   p is dim-by-np. Intended for dim >= 2.

if ~isstruct(param)
    error("bnd_cylinder: param must be a struct.");
end
if ~isfield(param, "center") || ~isfield(param, "axis") || ~isfield(param, "radius")
    error("bnd_cylinder: param must include center, axis, and radius.");
end

dim = size(p,1);
tol = 1e-8;
mode = "boundary";
if isfield(param, "tol")
    tol = param.tol;
end
if isfield(param, "mode")
    mode = param.mode;
end
if nargin >= 3
    mode = varargin{1};
end

center = param.center(:);
axisv = param.axis(:);
radius = param.radius;
finite = false;
tmin = -inf;
tmax = inf;

if numel(center) ~= dim || numel(axisv) ~= dim
    error("bnd_cylinder: center and axis length must equal size(p,1).");
end
if radius <= 0
    error("bnd_cylinder: radius must be positive.");
end

if isfield(param, "tmin") || isfield(param, "tmax")
    if ~(isfield(param, "tmin") && isfield(param, "tmax"))
        error("bnd_cylinder: provide both tmin and tmax.");
    end
    finite = true;
    tmin = param.tmin;
    tmax = param.tmax;
elseif isfield(param, "length")
    finite = true;
    L = param.length;
    tmin = -0.5 * L;
    tmax = 0.5 * L;
elseif isfield(param, "height")
    finite = true;
    L = param.height;
    tmin = -0.5 * L;
    tmax = 0.5 * L;
end

if finite
    if ~(isnumeric(tmin) && isnumeric(tmax) && isscalar(tmin) && isscalar(tmax))
        error("bnd_cylinder: tmin/tmax (or derived limits) must be numeric scalars.");
    end
    if tmax <= tmin
        error("bnd_cylinder: require tmax > tmin.");
    end
end

na = norm(axisv);
if na == 0
    error("bnd_cylinder: axis must be non-zero.");
end
axisu = axisv / na;

q = p - center;
ta = axisu.' * q;
proj = axisu * ta;
radial = q - proj;
r = sqrt(sum(radial.^2, 1));

if finite
    axialLoose = (ta >= tmin - tol) & (ta <= tmax + tol);
    axialStrict = (ta > tmin + tol) & (ta < tmax - tol);
else
    axialLoose = true(1, size(p,2));
    axialStrict = true(1, size(p,2));
end

radialLoose = r <= radius + tol;
radialStrict = r < radius - tol;
sideBoundary = axialLoose & (abs(r - radius) <= tol);
if finite
    capBoundary = radialLoose & ((abs(ta - tmin) <= tol) | (abs(ta - tmax) <= tol));
    onBoundary = sideBoundary | capBoundary;
else
    onBoundary = sideBoundary;
end

insideLoose = axialLoose & radialLoose;
insideStrict = axialStrict & radialStrict;
outsideLoose = (~axialStrict) | (r >= radius - tol);
outsideStrict = (~axialLoose) | (r > radius + tol);

m = applyMode(onBoundary, insideLoose, insideStrict, outsideLoose, outsideStrict, mode);
end

function m = applyMode(onBoundary, insideLoose, insideStrict, outsideLoose, outsideStrict, mode)
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
    error("bnd_cylinder: mode must be numeric, char, or string.");
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
    error("bnd_cylinder: unsupported mode '%s'.", ms);
end
end
