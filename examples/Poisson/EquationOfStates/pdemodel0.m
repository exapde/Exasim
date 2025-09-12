
function pde = pdemodel
pde.flux = @flux;
pde.source = @source;
pde.fbouhdg = @fbouhdg;
pde.initu = @initu;
end

function f = flux(u, q, w, v, x, t, mu, eta)
kappa = mu(1)*(1 + log(1 + u*u));
f = kappa*q;
end

function s = source(u, q, w, v, x, t, mu, eta)
x1 = x(1);
x2 = x(2);
s = x1*sin(x2);
end

function u0 = initu(x, mu, eta)
u0 = sym(1.0);
end

function fb = fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau)
fb1 = tau*(mu(2) - uhat);
f = flux(uhat, q, w, v, x, t, mu, eta);
fb2 = f(1)*n(1) + f(2)*n(2) + tau*(u(1)-uhat(1));
fb = [fb1 fb2];
end


