function mesh = mkmesh_platefluid(porder)

if nargin < 1
    porder = 2;
end

plateLength = 1.0;
plateThickness = 0.05;

plateX0 = 0.0;
plateX1 = plateLength;
plateY0 = 0.0;
plateY1 = plateThickness;

fluidWidth = 10.0;
fluidHeight = 10.0;

fluidXmin = 0.5 * (plateX0 + plateX1) - 0.5 * fluidWidth;
fluidXmax = fluidXmin + fluidWidth;
fluidYmin = 0.5 * (plateY0 + plateY1) - 0.5 * fluidHeight;
fluidYmax = fluidYmin + fluidHeight;

outer = [fluidXmin fluidYmin;
         fluidXmax fluidYmin;
         fluidXmax fluidYmax;
         fluidXmin fluidYmax];

inner = [plateX0 plateY0;
         plateX1 plateY0;
         plateX1 plateY1;
         plateX0 plateY1];

[p,t] = polymesh({outer, inner}, [1, 1], [1, 0; 0, 1], [0.60, 1.35], @href);
[p,t] = fixmesh(p,t);
elemtype = 0;
bndexpr = {'true'};
mesh = mkmesh(p, t, porder, bndexpr, elemtype, 1);

mesh.p = mesh.p';
mesh.t = mesh.t';

tol = 1e-6;
mesh.boundaryexpr = {...
    @(p) abs(p(2,:) - fluidYmin) < tol, ...
    @(p) abs(p(1,:) - fluidXmax) < tol, ...
    @(p) abs(p(2,:) - fluidYmax) < tol, ...
    @(p) abs(p(1,:) - fluidXmin) < tol, ...
    @(p) (p(1,:) >= plateX0 - tol) & (p(1,:) <= plateX1 + tol) & ...
         (p(2,:) >= plateY0 - tol) & (p(2,:) <= plateY1 + tol)};
mesh.boundarycondition = [2; 2; 2; 2; 2];
mesh.curvedboundary = [0 0 0 0 0];
mesh.curvedboundaryexpr = {...
    @(p) 0 * p(1,:), @(p) 0 * p(1,:), @(p) 0 * p(1,:), ...
    @(p) 0 * p(1,:), @(p) 0 * p(1,:)};
mesh.periodicexpr = {};
mesh.f = facenumbering(mesh.p, mesh.t, mesh.elemtype, mesh.boundaryexpr, mesh.periodicexpr);
mesh.xpe = mesh.plocal;
mesh.telem = mesh.tlocal;

figure(1); clf; meshplot(mesh);

function h = href(p)
dx = max(max(plateX0 - p(:,1), 0.0), p(:,1) - plateX1);
dy = max(max(plateY0 - p(:,2), 0.0), p(:,2) - plateY1);
d = sqrt(dx.^2 + dy.^2);

h = 0.60 * ones(size(d));
h(d < 1.50) = 0.20;
h(d < 0.60) = 0.08;
h(d < 0.15) = 0.04;
end

end
