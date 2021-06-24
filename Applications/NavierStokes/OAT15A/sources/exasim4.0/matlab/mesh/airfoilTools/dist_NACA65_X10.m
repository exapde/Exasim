
function dist = dist_NACA65_X10(p,X,a)
% Distance function for NACA 65-010 airfoil.
% (piecewise linear approximation)

if nargin<3; a = 1.0; end

[xf,yf] = NACA65_X10(X,a);
pf = [xf yf];
dist = dpoly(p,pf);

end
