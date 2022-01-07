function [p,t] = lesmesh2d(xf, yf, dlay, dwall, nx, ny, xref, yref)

% (xf, yf): airfoil coordinates
% dlay: BL thickness
% dwall: thickness of the first element at the wall
% nx: number of elements in x-direction
% ny: number of elements in y-direction
% xref: refinement along the x-direction
% yref: refinement along the y-direction
  
% Example: 
% dlay=0.1; dwall=2e-5; nx=96; ny = 25;
% [p,t] = lesmesh2d(xf, yf, dlay, dwall, nx, ny, [1.3 2; 1.55 1.78], [0.03 0.01 0.003]);

if size(xref,2)~=2
    error('xref must have dimension Nx times 2');
end

% calculate the mesh ratio
c = 1 - dlay/dwall;
rat = fsolve(@(x) scalingfun(x,ny,c),[1;3]);
rat = rat(1); 

% scaling distribution over the normal direction
yv = zeros(ny+1,1);
yv(2) = dwall;
for i = 1:(ny-1)
    yv(i+2) = yv(i+1) + dwall*(rat^i);
end

if abs(yv(end)-dlay)>1e-8
    error("Something wrong with the input parameters (dlay, dwall, ny)");
end

% Uniform distribution over foil
ns = length(xf);
t = 0:ns-1;
spx = spline(t,xf);
spy = spline(t,yf);
ttp = distribute(nx,spx,spy,ns);
xv = zeros(nx+1,1);
xv(1:nx+1,1) = 2*ttp/(ns-1);

% refine according to xref
for i = 1:size(xref,1)
    ind1 = (xv < xref(i,1));
    ind2 = (xv >= xref(i,1)) & (xv <= xref(i,2));
    ind3 = (xv > xref(i,2));
    x1 = xv(ind1); % no refinement
    x2 = xv(ind2); % yes refinement
    x3 = xv(ind3); % no refinement
    % refine x2
    xw = 0.5*(x2(1:end-1)+x2(2:end));
    x2 = sort([x2; xw]);
    % merge grid points 
    xv = unique([x1; x2; x3]);
end

size(xv)

% make the grid from points
[p,t] = quadgrid(xv,yv);

size(t)

% refine according to yref
n = length(yref);
if n>0
    yref = sort(yref,'descend');
    for i = 1:n
        [p,t] = refineaty(p,t,yref(i));
    end
    [p,t] = fixmesh(p,t);
end

[p,t] = removeelemement(p, t, ['y>' num2str(dlay/1.75)]);

% plot grid 
figure(1);clf;simpplot(p,t);axis on;

% map rect to foil
p = map2foil(p,xf,yf);

% plot grid 
figure(2);clf;simpplot(p,t);axis on;


