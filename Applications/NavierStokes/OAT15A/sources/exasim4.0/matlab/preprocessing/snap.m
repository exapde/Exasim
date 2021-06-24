function x=snap(x)

tol=sqrt(eps);
x=tol*round(x/tol);
