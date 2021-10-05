function [p,t] = lesmesh3d_box(dlay, dwall, ny, xv, zv, yref)

% xv: grid points along z-direction
% zv: grid points along z-direction
% dlay: BL thickness
% dwall: thickness of the first element at the wall
% ny: number of elements in y-direction
% yref: refinement along the y-direction
  
% % Example: 
% dlay=0.06; dwall=2e-5; nx=48; ny = 30; zv = linspace(0,0.2,11);
% [p,t] = lesmesh3d(xf, yf, zv, dlay, dwall, nx, ny, [0.6 2; 0.9 2; 1.3 1.85], [0.045 0.025 0.01]);
% porder = 1; bndexpr = {'true'}; elemtype = 1; nodetype=0;
% mesh = mkmesh(p,t,porder,bndexpr,elemtype,nodetype);

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

% make the grid from points
[p,t] = hexgrid(xv,yv,zv);

% refine according to yref
yref = sort(yref,'descend');
for i = 1:length(yref)
    [p,t] = refineaty3d2to1(p,t,yref(i));
end


