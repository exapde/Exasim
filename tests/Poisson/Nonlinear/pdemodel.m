
function pde = pdemodel
pde.flux = @flux;
pde.source = @source;
pde.fbou = @fbou;
pde.ubou = @ubou;
pde.initu = @initu;
end

function f = flux(u, q, w, v, x, t, mu, eta)
f = (1+u*u)*mu(1)*q;
end

function s = source(u, q, w, v, x, t, mu, eta)
x1 = x(1);
x2 = x(2);
s = x1*sin(x2);
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
f = flux(u, q, w, v, x, t, mu, eta);
tm = f(1)*n(1) + f(2)*n(2) + tau*(u(1)-uhat(1));
fb = [tm 0.0]; % [Dirichlet, Neumman]
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
ub = [mu(2) u]; % [Dirichlet, Neumman]
end

function u0 = initu(x, mu, eta)
u0 = sym(1.0);
end

