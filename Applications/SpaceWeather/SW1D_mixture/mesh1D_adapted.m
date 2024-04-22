function [p,t] = mesh1D_adapted(r1,r2,ndiv)

dlay = 1;
dwall = 5e-3; 

% calculate the mesh ratio
c = 1 - dlay/dwall;
rat = fsolve(@(x) scalingfun(x,ndiv,c),[1;1.5]);
rat = rat(1);

% scaling distribution over the normal direction
xv = zeros(ndiv+1,1);
xv(2) = dwall;
for i = 1:(ndiv-1)
    xv(i+2) = xv(i+1) + dwall*(rat^i);
end

if abs(xv(end)-dlay)>1e-8
    error("Something wrong with the input parameters (dlay, dwall, nx)");
end

p = xv'*(r2-r1) + r1;
t = [(1:ndiv); (2:ndiv+1)];
