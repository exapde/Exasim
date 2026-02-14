function m = bnd_sphere(p, param, varargin)
%BND_SPHERE Predicate for a sphere (or circle in 2D).
%   m = bnd_sphere(p, param) returns logical row vector m for points on:
%       ||p - center|| = radius (within tolerance).
%
%   m = bnd_sphere(p, param, mode) also supports logical selections:
%       mode = "boundary" (default), "inside", "outside",
%              "strict-inside", "strict-outside"
%
%   Supported parameterizations:
%   1) Numeric: [center(1:dim), radius] or [center(1:dim), radius, tol]
%   2) Struct:  param.center, param.radius, optional param.tol, optional param.mode

% % On sphere (default)
% @(p) bnd_sphere(p, [0, 0, 1.0, 1e-6])
% 
% % Inside sphere
% @(p) bnd_sphere(p, [0, 0, 1.0, 1e-6], "inside")
% 
% % Outside sphere
% @(p) bnd_sphere(p, struct("center",[0;0], "radius",1.0, "tol",1e-6, "mode","outside"))

dim = size(p,1);
tol = 1e-8;
mode = "boundary";

if isnumeric(param)
    if numel(param) < dim + 1
        error("bnd_sphere: numeric param must be [center(1:dim), radius, (tol)].");
    end
    center = param(1:dim).';
    radius = param(dim+1);
    if numel(param) >= dim + 2
        tol = param(dim+2);
    end
elseif isstruct(param)
    if ~isfield(param, "center") || ~isfield(param, "radius")
        error("bnd_sphere: struct param requires center and radius.");
    end
    center = param.center(:);
    radius = param.radius;
    if isfield(param, "tol")
        tol = param.tol;
    end
    if isfield(param, "mode")
        mode = param.mode;
    end
else
    error("bnd_sphere: param must be numeric or struct.");
end

if nargin >= 3
    mode = varargin{1};
end

if numel(center) ~= dim
    error("bnd_sphere: center length must equal size(p,1).");
end

r = sqrt(sum((p - center).^2, 1));
m = applyMode(r, radius, tol, mode);
end

function m = applyMode(r, radius, tol, mode)
if isnumeric(mode)
    s = sign(mode);
    if s == 0
        m = abs(r - radius) <= tol;
    elseif s > 0
        m = r <= radius + tol;
    else
        m = r >= radius - tol;
    end
    return;
end

if ~(ischar(mode) || isstring(mode))
    error("bnd_sphere: mode must be numeric, char, or string.");
end

ms = lower(strtrim(char(mode)));
if strcmp(ms, "boundary") || strcmp(ms, "on") || strcmp(ms, "surface") || strcmp(ms, "=")
    m = abs(r - radius) <= tol;
elseif strcmp(ms, "inside") || strcmp(ms, "in") || strcmp(ms, "interior") || strcmp(ms, "<") || strcmp(ms, "<=")
    m = r <= radius + tol;
elseif strcmp(ms, "outside") || strcmp(ms, "out") || strcmp(ms, "exterior") || strcmp(ms, ">") || strcmp(ms, ">=")
    m = r >= radius - tol;
elseif strcmp(ms, "strict-inside") || strcmp(ms, "inside-strict")
    m = r < radius - tol;
elseif strcmp(ms, "strict-outside") || strcmp(ms, "outside-strict")
    m = r > radius + tol;
else
    error("bnd_sphere: unsupported mode '%s'.", ms);
end
end
