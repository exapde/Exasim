function pde = pdemodel
pde.mass = @mass;
pde.flux = @flux;
pde.source = @source;
pde.fbou = @fbou;
pde.ubou = @ubou;
pde.initu = @initu;
pde.initq = @initq;
pde.initw = @initw;
end

function m = mass(u, q, w, v, x, t, mu, eta)
m = sym(1.0);
end

function f = flux(u, q, w, v, x, t, mu, eta)
f = mu*q;
end

function s = source(u, q, w, v, x, t, mu, eta)
s = sym(0.0);
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
f = flux(u, q, w, v, x, t, mu, eta);
fb = f(1)*n(1) + f(2)*n(2) + tau*(u(1)-0.0);
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
ub = sym(0.0); 
end

function u0 = initu(x, mu, eta)
x1 = x(1);
x2 = x(2);
u0 = sin(pi*x1)*sin(pi*x2);
end

function q0 = initq(x, mu, eta)
q0 = [sym(0.0) sym(0.0)];
end

function w0 = initw(x, mu, eta)
w0 = sym(0.0);
end


