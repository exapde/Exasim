function [p, t] = squaregrid(n, elemtype)

[p,t] = squaremesh(n,n,1,elemtype);
tm = [size(p) size(t) p(:)' t(:)']; writebin("grid.bin",tm);

