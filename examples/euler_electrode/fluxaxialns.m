function f = fluxaxialns(u, q, w, v, x, t, mu, eta)
   % pde.physicsmu = [gam Re Pr Minf rinf ruinf rvinf rEinf Tref avk avs];

% regularization mueters
alpha = 4.0e3;
rmin = 5.0e-2;
pmin = 2.0e-3;

r = u(1);
ru = u(2);
rv = u(3);
rE = u(4);
rx = q(1);
rux = q(2);
rvx = q(3);
rEx = q(4);
ry = q(5);
ruy = q(6);
rvy = q(7);
rEy = q(8);

% Regularization of rho (cannot be smaller than rmin)
r = rmin + lmax(r-rmin,alpha);
% Density sensor
dr = atan(alpha*(r - rmin))/pi + (alpha*(r - rmin))/(pi*(alpha^2*(r - rmin)^2 + 1)) + 1/2;
r1 = 1/r;
uv = ru*r1;
vv = rv*r1;
E = rE*r1;

gam = 1.4;        % Thermodynamic quantities are determined consistently from the equation of state and used to couple the energy to the momentum equation in the following manner: energy->temperature->pressure
e = E-0.5*(uv*uv+vv*vv);
T = e*(gam-1);
p = r*T;
H = E + T;

p = pmin + lmax(p-pmin,alpha);

% Inviscid fluxes only -- Euler
fi = [ru, ru*uv+p, rv*uv, ru*H, ...
      rv, ru*vv, rv*vv+p, rv*H];

f = fi;

f = reshape(f,[4,2]);        

end

% % Perform a nonlinear solve for the temperature given the Cp function
% function cp = cpMassAir(T)

%       if and(T >= 200, T < 1000)
%       cp = 8.06894772E+04*T^-2 -1.58978421E+03*T^-1 + 4.14953775E+01 + -4.74713424E-02*T + 8.85482153E-05*T^2 -6.61274683E-08*T^3 + 1.82416968E-11*T^4;
      
%       elseif and(T >= 1000, T < 6000)
%       cp = 1.91291177E+06*T^-2 -1.02330966E+04*T^-1 + 4.26747334E+01 -1.66264426E-03*T + 5.67244351E-07*T^2 -8.71530239E-11*T^3 + 5.38836908E-15*T^4;
      
%       elseif and(T >= 6000, T <= 20000)
%       cp = 6.29946629E+09*T^-2 -4.68831578E+06*T^-1 + 1.43258522E+03 -2.10067001E-01*T + 1.66801152E-05*T^2 -6.43037167E-10*T^3 + 9.45602828E-15*T^4;
%       else
%       end
% end
