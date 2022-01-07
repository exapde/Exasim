function f = limiting(x,xmin,xmax,alpha,beta)

% f = xmin + max(x-beta,0);
% f = min(f-xmax,0)+xmax;

f = xmin + lmax(x-beta,alpha);
f = lmin(f-xmax,alpha) + xmax;

