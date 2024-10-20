
function pde = pdemodel
pde.mass = @mass;
pde.flux = @flux;
pde.source = @source;
pde.fbou = @fbou;
pde.fbouhdg = @fbouhdg;
pde.ubou = @ubou;
pde.initu = @initu;
end

function m = mass(u, q, w, v, x, t, mu, eta)
m = sym(1.0);
end

function f = flux(u, q, w, v, x, t, mu, eta)
f = [mu(1)*u(1) mu(2)*u(1)];
end

function s = source(u, q, w, v, x, t, mu, eta)
s = sym(0.0);
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
f = flux(u, q, w, v, x, t, mu, eta);
fb = f(1)*n(1) + f(2)*n(2) + tau*(u(1)-uhat(1));
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
ub = sym(0.0); 
end

function u0 = initu(x, mu, eta)
u0 = exp(-(x(1)-0.5)*(x(1)-0.5)/0.02-(x(2)-0.5)*(x(2)-0.5)/0.02);
end

function fb = fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau)
fb = tau*(sym(0.0)-uhat(1));
end
