dlay = 1;
dwall = 3e-3; 
nx = 240;
ny = 24; 
xv = linspace(0, 1, nx); 
yref = [3e-2 12e-2 3e-1 6e-1];

[p,t,yv] = lesmesh2d_rect(dlay, dwall, ny, xv, yref);

p = p';
r1 = 1;
r2 = r1 + 0.1;
pnew = p;
pnew(:,1) = -(r1+(r2-r1)*p(:,2)).*sin(2*pi*p(:,1));
pnew(:,2) = -(r1+(r2-r1)*p(:,2)).*cos(2*pi*p(:,1));
p = pnew';

simpplot(p,t); axis equal; axis tight;
