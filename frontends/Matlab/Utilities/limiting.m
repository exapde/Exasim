function f = limiting(x,fmin,fmax,alpha,beta)

% f = fmin + max(x-beta,0);
% f = min(f-fmax,0)+fmax;

f = fmin + lmax(x-beta,alpha);
f = lmin(f-fmax,alpha) + fmax;

