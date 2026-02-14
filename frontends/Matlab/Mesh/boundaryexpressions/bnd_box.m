function m = bnd_box(p, param, varargin)
%BND_BOX Predicate for an axis-aligned box.
%   m = bnd_box(p, param) returns logical row vector m for points on any
%   box face (default mode), within tolerance.
%
%   m = bnd_box(p, param, mode) supports:
%       mode = "boundary" (default), "inside", "outside",
%              "strict-inside", "strict-outside"
%
%   Supported parameterizations:
%   1) Numeric: [xmin xmax ymin ymax ...] or [xmin xmax ymin ymax ... tol]
%      (2*dim values, optional trailing tol)
%   2) Struct:  param.min, param.max, optional param.tol, optional param.mode
%
%   p is dim-by-np.

dim = size(p,1);
tol = 1e-8;
mode = "boundary";

if isnumeric(param)
    if numel(param) ~= 2*dim && numel(param) ~= 2*dim + 1
        error("bnd_box: numeric param must have 2*dim values (+ optional tol).");
    end
    if numel(param) == 2*dim + 1
        tol = param(end);
        vals = param(1:end-1);
    else
        vals = param;
    end
    b = reshape(vals, 2, []).';
    bmin = b(:,1);
    bmax = b(:,2);
elseif isstruct(param)
    if ~isfield(param, "min") || ~isfield(param, "max")
        error("bnd_box: struct param requires min and max.");
    end
    bmin = param.min(:);
    bmax = param.max(:);
    if isfield(param, "tol")
        tol = param.tol;
    end
    if isfield(param, "mode")
        mode = param.mode;
    end
else
    error("bnd_box: param must be numeric or struct.");
end

if nargin >= 3
    mode = varargin{1};
end

if numel(bmin) ~= dim || numel(bmax) ~= dim
    error("bnd_box: min/max length must equal size(p,1).");
end

insideLoose = all((p >= (bmin - tol)) & (p <= (bmax + tol)), 1);
insideStrict = all((p > (bmin + tol)) & (p < (bmax - tol)), 1);
outsideLoose = any((p <= (bmin + tol)) | (p >= (bmax - tol)), 1);
outsideStrict = any((p < (bmin - tol)) | (p > (bmax + tol)), 1);
onface = false(1, size(p,2));

for i = 1:dim
    onface = onface | (abs(p(i,:) - bmin(i)) <= tol) | (abs(p(i,:) - bmax(i)) <= tol);
end

m = selectMode(mode, insideLoose, insideStrict, outsideLoose, outsideStrict, onface);
end

function m = selectMode(mode, insideLoose, insideStrict, outsideLoose, outsideStrict, onface)
if isnumeric(mode)
    s = sign(mode);
    if s == 0
        m = insideLoose & onface;
    elseif s > 0
        m = insideLoose;
    else
        m = outsideLoose;
    end
    return;
end

if ~(ischar(mode) || isstring(mode))
    error("bnd_box: mode must be numeric, char, or string.");
end

ms = lower(strtrim(char(mode)));
if strcmp(ms, "boundary") || strcmp(ms, "on") || strcmp(ms, "surface") || strcmp(ms, "=")
    m = insideLoose & onface;
elseif strcmp(ms, "inside") || strcmp(ms, "in") || strcmp(ms, "interior") || strcmp(ms, "<") || strcmp(ms, "<=")
    m = insideLoose;
elseif strcmp(ms, "outside") || strcmp(ms, "out") || strcmp(ms, "exterior") || strcmp(ms, ">") || strcmp(ms, ">=")
    m = outsideLoose;
elseif strcmp(ms, "strict-inside") || strcmp(ms, "inside-strict")
    m = insideStrict;
elseif strcmp(ms, "strict-outside") || strcmp(ms, "outside-strict")
    m = outsideStrict;
else
    error("bnd_box: unsupported mode '%s'.", ms);
end
end
