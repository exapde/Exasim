function [p,t,yv] = lesmesh2d_rect(dlay, dwall, ny, xv, yref)
% Generate an adaptive boundary layer mesh for the flat plate

% dlay:  boundary layer thickness (in y-direction)
% dwall: thickness of the first element at the wall
% ny:    number of elements in y-direction
% xv:    points along x-direction
% yref:  points along y-direction at which the mesh is refined
  
% Example: 
% dlay=0.1; dwall=2e-5; ny = 25; xv = linspace(0, 1, 1000); yref = [1e-4 1e-3 1e-2];
% [p,t,yv] = lesmesh2d_rect(dlay, dwall, ny, xv, yref);

% calculate the mesh ratio
c = 1 - dlay/dwall;
rat = fsolve(@(x) scalingfun(x,ny,c),[1;1.5]);
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

% make the grid from points
[p,t] = quadgrid(xv,yv);

% refine according to yref
n = length(yref);
if n>0
    yref = sort(yref,'descend');
    for i = 1:n
        [p,t] = refineaty(p,t,yref(i));
    end
    [p,t] = fixmesh(p,t);
end

p = p';
t = t';

% plot grid 
%figure(1);clf;simpplot(p,t);axis on; axis equal;




