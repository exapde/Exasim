function f = lmax(x,alpha)

f = x.*(atan(alpha*x)/pi + 0.5) - atan(alpha)/pi + 0.5;
