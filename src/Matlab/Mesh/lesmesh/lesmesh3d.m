function [p,t] = lesmesh3d(xf, yf, zv, dlay, dwall,  nx, ny, xref, yref)

% (xf, yf): airfoil coordinates
% zv: grid points along z-direction
% dlay: BL thickness
% dwall: thickness of the first element at the wall
% nx: number of elements in x-direction
% ny: number of elements in y-direction
% xref: refinement along the x-direction
% yref: refinement along the y-direction
  
% Example: 
% dlay=0.06; dwall=2e-5; nx=48; ny = 30; zv = linspace(0,0.2,11);
% [p,t] = lesmesh3d(xf, yf, zv, dlay, dwall, nx, ny, [0.6 2; 0.9 2; 1.3 1.85], [0.045 0.025 0.01]);
% porder = 1; bndexpr = {'true'}; elemtype = 1; nodetype=0;
% mesh = mkmesh(p,t,porder,bndexpr,elemtype,nodetype);

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

% make the grid from points
[p,t] = hexgrid(xv,yv,zv);

% refine according to yref
yref = sort(yref,'descend');
for i = 1:length(yref)
    [p,t] = refineaty3d(p,t,yref(i));
end

p(:,1:2) = map2foil(p(:,1:2),xf,yf);

% % refine according to yref
% yref = sort(yref,'descend');
% for i = 1:length(yref)
%     [p,t] = refineaty(p,t,yref(i));
% end
% [p,t] = fixmesh(p,t);
% 
% % plot grid 
% figure(1);clf;simpplot(p,t);axis on;
% 
% % map rect to foil
% p = map2foil(p,xf,yf);
% 
% % plot grid 
% figure(2);clf;simpplot(p,t);axis on;
% 

