function [uo] = get_outletBC(u,param,nlg)


% Get the regularization parameters
rmin = param(5);
pmin = param(7);
gam  = param(1);
alpha = 1.0e3;

% Get the Target pressure
pout = param(11);

% Outwards normals
nx = nlg(1); 
ny = nlg(2); 
nz = nlg(3);

% Regularization of rho and pressure
r = rmin + lmax(u(1)-rmin,alpha);
uv = u(2)/r; 
vv = u(3)/r;
wv = u(3)/r;
p = (gam-1)*(u(4) - 0.5*(u(2)*uv + u(3)*vv + u(4)*wv));
p = pmin + lmax(p-pmin,alpha);

% Compute the local Mach number in the normal direction
vn = lmax((uv*nx + vv*ny + wv*nz),alpha);
Ma = vn / sqrt(gam*p/r);
% step-ification of Ma :
alpha = 100;
Mastep = 0.5*(1+tanh(alpha*(Ma-1)));

% Subsonic outlet variables
usub = u;
usub(1) = u(1); % + (u(1)/p)*(pout-p); %u(1); %pout/p * u(1);
usub(2) = usub(1)*uv;
usub(3) = usub(1)*vv;
usub(4) = usub(1)*wv;
usub(5) = pout/(gam-1) + 0.5*usub(2)*usub(2)/usub(1) + 0.5*usub(3)*usub(3)/usub(1) + 0.5*usub(4)*usub(4)/usub(1);

% Subsonic outlet BC
uo = usub*(1-Mastep) + u*Mastep;


end

