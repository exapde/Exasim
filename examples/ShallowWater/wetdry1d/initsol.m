function [H0, Hu0] = initsol(x, mu)
hfloor = mu(7);
[zb, ~] = bathymetry(x, mu);
H0 = positivepart(-zb);
Hu0 = 0*x;
end

function H = positivepart(h)
gamma = 1e6;
H = h.*(atan(gamma*(h))/pi + 0.5) - atan(gamma)/pi + 0.5;    
H = H + min(H(:));
end

function [zb, dzbdx] = bathymetry(x, mu)
L = mu(2);
zbeach = mu(3);
zdeep = mu(4);

slope = (zdeep - zbeach)/L;
zb = zbeach + slope*x;
dzbdx = slope;
end

