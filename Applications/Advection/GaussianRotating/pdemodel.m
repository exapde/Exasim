% Poisson equation with homogenous Dirichlet condition on a unit square
% q + \nabla u = 0 in \Omega
% \nabla dot flux = source in \Omega
% u = 0 on \partial Omega
% flux = 2*param*q
% source = 2*pi*pi*sin(pi*x1)*sin(pi*x2);

function pde = pdemodel
pde.mass = @mass;
pde.flux = @flux;
pde.source = @source;
pde.fbou = @fbou;
pde.ubou = @ubou;
pde.initu = @initu;
end

function m = mass(u, q, w, v, x, t, mu, eta)
m = sym(1.0);
end

function f = flux(u, q, w, v, x, t, mu, eta)
f = [mu(1)*x(2)*u(1) -mu(2)*x(1)*u(1)];
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
u0 = exp(-(x(1)-0.25)*(x(1)-0.25)/0.01-(x(2)-0.0)*(x(2)-0.0)/0.01);
end

