
function pde = pdemodel
pde.flux = @flux;
pde.source = @source;
pde.sourcew = @sourcew;
pde.fbou = @fbou;
pde.fbouhdg = @fbouhdg;
pde.ubou = @ubou;
pde.initu = @initu;
pde.initw = @initw;
pde.eos = @eos;
end

function f = flux(u, q, w, v, x, t, mu, eta)
kappa = mu(1)*(1 + w);
f = kappa*q;
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

function w0 = initw(x, mu, eta)
w0 = sym(2.3026);
end

function s = eos(u, q, w, v, x, t, mu, eta)
s = exp(w) - sym(1.0) - u*u; 
end

function fb = fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau)
fb1 = tau*(mu(2) - uhat);
f = flux(uhat, q, w, v, x, t, mu, eta);
fb2 = f(1)*n(1) + f(2)*n(2) + tau*(u(1)-uhat(1));
fb = [fb1 fb2];
end

function s = sourcew(u, q, w, v, x, t, mu, eta)
s = exp(w) - sym(1.0) - u*u; 
end


