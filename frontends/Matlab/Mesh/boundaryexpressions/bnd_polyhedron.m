function m = bnd_polyhedron(p, param, varargin)
%BND_POLYHEDRON Predicate for a convex polyhedron (halfspace intersection).
%   m = bnd_polyhedron(p, param) returns points on the polyhedron boundary.
%
%   m = bnd_polyhedron(p, param, mode) supports:
%       mode = "boundary" (default), "inside", "outside",
%              "strict-inside", "strict-outside"
%
%   Supported parameterizations:
%   1) A*x <= b form:
%      param.A : nf-by-dim matrix
%      param.b : nf-by-1 (or 1-by-nf) vector
%   2) normals/offsets form:
%      param.normals : dim-by-nf or nf-by-dim
%      param.offsets : nf-by-1 (or 1-by-nf)
%      Interpreted as normals(:,i).' * x + offsets(i) <= 0.
%
%   Optional struct fields: param.tol (default 1e-8), param.mode.
%   p is dim-by-np.

if ~isstruct(param)
    error("bnd_polyhedron: param must be a struct.");
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

[A, b] = parseConstraints(param, dim);

g = A * p - b; % inside if g <= 0

insideLoose = all(g <= tol, 1);
insideStrict = all(g < -tol, 1);
onBoundary = insideLoose & any(abs(g) <= tol, 1);
outsideLoose = ~insideStrict;
outsideStrict = any(g > tol, 1);

m = selectMode(onBoundary, insideLoose, insideStrict, outsideLoose, outsideStrict, mode);
end

function [A, b] = parseConstraints(param, dim)
if isfield(param, "A") && isfield(param, "b")
    A = param.A;
    b = param.b(:);
    if ~isnumeric(A) || ~isnumeric(b)
        error("bnd_polyhedron: A and b must be numeric.");
    end
    if size(A,2) ~= dim
        error("bnd_polyhedron: A must be nf-by-size(p,1).");
    end
    if size(A,1) ~= numel(b)
        error("bnd_polyhedron: length(b) must match size(A,1).");
    end
    return;
end

if isfield(param, "normals") && isfield(param, "offsets")
    N = param.normals;
    c = param.offsets(:);
    if ~isnumeric(N) || ~isnumeric(c)
        error("bnd_polyhedron: normals and offsets must be numeric.");
    end

    if size(N,1) == dim
        if size(N,2) ~= numel(c)
            error("bnd_polyhedron: number of normals must match offsets.");
        end
        A = N.';
    elseif size(N,2) == dim
        if size(N,1) ~= numel(c)
            error("bnd_polyhedron: number of normals must match offsets.");
        end
        A = N;
    else
        error("bnd_polyhedron: normals must be dim-by-nf or nf-by-dim.");
    end

    b = -c;
    return;
end

error("bnd_polyhedron: provide either (A,b) or (normals,offsets).");
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
    error("bnd_polyhedron: mode must be numeric, char, or string.");
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
    error("bnd_polyhedron: unsupported mode '%s'.", ms);
end
end
