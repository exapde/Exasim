function [p,t,yv] = lesmesh2d_rect(dlay, dwall, ny, xv, yref)

% dlay: BL thickness
% dwall: thickness of the first element at the wall
% ny: number of elements in y-direction
% xv: points along x-direction
  
% Example: 
% dlay=0.1; dwall=2e-5; nx=96; ny = 25;
% [p,t] = lesmesh2d(xf, yf, dlay, dwall, nx, ny, [1.3 2; 1.55 1.78], [0.03 0.01 0.003]);

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

p=p';
t=t';

% plot grid 
% figure(1);clf;simpplot(p,t);axis on;



