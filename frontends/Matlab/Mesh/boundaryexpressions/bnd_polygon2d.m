function m = bnd_polygon2d(p, param, varargin)
%BND_POLYGON2D Predicate for a 2D polygon region.
%   m = bnd_polygon2d(p, param) returns points on polygon boundary (default).
%
%   m = bnd_polygon2d(p, param, mode) supports:
%       mode = "boundary" (default), "inside", "outside",
%              "strict-inside", "strict-outside"
%
%   Required fields / data:
%   - Numeric param: vertices as nv-by-2 or 2-by-nv.
%   - Struct param:  param.vertices (same formats).
%
%   Optional struct fields:
%   - param.dims : two coordinate indices from p (default [1 2])
%   - param.tol  : boundary tolerance (default 1e-8)
%   - param.mode : selection mode
%
%   p is dim-by-np.

tol = 1e-8;
mode = "boundary";
dims = [1 2];

if isnumeric(param)
    verts = parseVertices(param);
elseif isstruct(param)
    if ~isfield(param, "vertices")
        error("bnd_polygon2d: struct param requires vertices.");
    end
    verts = parseVertices(param.vertices);
    if isfield(param, "dims")
        dims = param.dims(:).';
    end
    if isfield(param, "tol")
        tol = param.tol;
    end
    if isfield(param, "mode")
        mode = param.mode;
    end
else
    error("bnd_polygon2d: param must be numeric or struct.");
end

if nargin >= 3
    mode = varargin{1};
end

if numel(dims) ~= 2
    error("bnd_polygon2d: dims must contain two indices.");
end
dims = round(dims);
if any(dims < 1) || any(dims > size(p,1))
    error("bnd_polygon2d: dims must be between 1 and size(p,1).");
end

xy = p(dims, :);
x = xy(1,:);
y = xy(2,:);
vx = verts(1,:);
vy = verts(2,:);

[in0, on0] = inpolygon(x, y, vx, vy);
dist = pointPolygonDistance(xy, verts);
onBoundary = on0 | (dist <= tol);
insideLoose = in0 | (dist <= tol);
insideStrict = in0 & (~onBoundary);
outsideLoose = ~insideStrict;
outsideStrict = (~in0) & (dist > tol);

m = selectMode(onBoundary, insideLoose, insideStrict, outsideLoose, outsideStrict, mode);
end

function verts = parseVertices(v)
if isempty(v) || ~isnumeric(v)
    error("bnd_polygon2d: vertices must be a numeric array.");
end

if size(v,2) == 2
    verts = v.';
elseif size(v,1) == 2
    verts = v;
else
    error("bnd_polygon2d: vertices must be nv-by-2 or 2-by-nv.");
end

if size(verts,2) < 3
    error("bnd_polygon2d: polygon requires at least 3 vertices.");
end
end

function dist = pointPolygonDistance(xy, verts)
np = size(xy,2);
nv = size(verts,2);
dist2 = inf(1,np);

for i = 1:nv
    i2 = i + 1;
    if i == nv
        i2 = 1;
    end
    a = verts(:,i);
    b = verts(:,i2);
    d2 = pointSegmentDistance2(xy, a, b);
    dist2 = min(dist2, d2);
end

dist = sqrt(dist2);
end

function d2 = pointSegmentDistance2(xy, a, b)
ab = b - a;
ab2 = sum(ab.^2);
if ab2 == 0
    d = xy - a;
    d2 = sum(d.^2, 1);
    return;
end

t = (ab.' * (xy - a)) / ab2;
t = max(0, min(1, t));
q = a + ab * t;
d = xy - q;
d2 = sum(d.^2, 1);
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
    error("bnd_polygon2d: mode must be numeric, char, or string.");
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
    error("bnd_polygon2d: unsupported mode '%s'.", ms);
end
end
