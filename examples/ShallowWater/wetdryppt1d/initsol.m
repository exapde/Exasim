function [phi0, H0] = initsol(x, mu)
Hmin = mu(7);
%Heps = mu(8);

zb = localbathymetry(x, mu);
%H0 = smoothmax(-zb, Hmin + Heps, Heps);
gamma = 1e2;
h = - zb;
H0 = h.*(atan(gamma*(h))/pi + 0.5) - atan(gamma)/pi + 0.5 + Hmin;    
phi0 = log(H0);
end

% function y = smoothmax(a, b, eps0)
% y = 0.5*(a + b + sqrt((a - b).*(a - b) + eps0*eps0));
% end
