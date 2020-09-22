% Poisson equation with homogenous Dirichlet condition on a unit square
% q + \nabla u = 0 in \Omega
% \nabla dot flux = source in \Omega
% u = 0 on \partial Omega
% flux = 2*param*q
% source = 2*pi*pi*sin(pi*x1)*sin(pi*x2);

function pde = pdemodel
pde.flux = @flux;
pde.source = @source;
pde.fbou = @fbou;
pde.ubou = @ubou;
pde.initu = @initu;
end

function f = flux(u, q, w, v, x, t, mu, eta)
f = mu(1)*q + u(1)*mu(2:3);
end

function s = source(u, q, w, v, x, t, mu, eta)
x1 = x(1);
x2 = x(2);
v1 = mu(2);
v2 = mu(3);
s = (v1*x2*exp(v1)*exp(v2)*(exp(v2*(x2 - 1)) - 1)*(exp(v1*(x1 - 1)) + v1*x1*exp(v1*(x1 - 1)) - 1))/((exp(v1) - 1)*(exp(v2) - 1)) - (v2*x1*exp(v1 + v2*x2)*(exp(v1*(x1 - 1)) - 1)*(v2*x2 + 2))/((exp(v1) - 1)*(exp(v2) - 1)) - (v1*x2*exp(v2 + v1*x1)*(exp(v2*(x2 - 1)) - 1)*(v1*x1 + 2))/((exp(v1) - 1)*(exp(v2) - 1)) + (v2*x1*exp(v1)*exp(v2)*(exp(v1*(x1 - 1)) - 1)*(exp(v2*(x2 - 1)) + v2*x2*exp(v2*(x2 - 1)) - 1))/((exp(v1) - 1)*(exp(v2) - 1)); 
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
f = flux(u, q, w, v, x, t, mu, eta);
fb = f(1)*n(1) + f(2)*n(2) + tau*(u(1)-uhat(1));
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
ub = sym(0.0); 
end

function u0 = initu(x, mu, eta)
u0 = sym(0.0);
end

