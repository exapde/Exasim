function f = evfunc(x, y, a, b)

f = 0*x;
ind = find(y <= a*x);
f(ind) = b * y(ind) .* (a*x(ind) - y(ind)).^2;


