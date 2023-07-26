function [p,t] = mesh1D_adapted(r1,r2)

dlay = 1;
dwall = 1e-2; 
nx = 16;

% calculate the mesh ratio
c = 1 - dlay/dwall;
rat = fsolve(@(x) scalingfun(x,nx,c),[1;1.5]);
rat = rat(1);

% scaling distribution over the normal direction
xv = zeros(nx+1,1);
xv(2) = dwall;
for i = 1:(nx-1)
    xv(i+2) = xv(i+1) + dwall*(rat^i);
end

if abs(xv(end)-dlay)>1e-8
    error("Something wrong with the input parameters (dlay, dwall, nx)");
end

p = xv'*(r2-r1) + r1;
t = [(1:nx); (2:nx+1)];
