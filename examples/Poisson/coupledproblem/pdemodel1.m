
function pde = pdemodel1
pde.flux = @flux;
pde.source = @source;
pde.fbou = @fbou;
pde.fint = @fint;
pde.fbouhdg = @fbouhdg;
pde.ubou = @ubou;
pde.initu = @initu;
end

function f = flux(u, q, w, v, x, t, mu, eta)
f = mu*q;
end

function s = source(u, q, w, v, x, t, mu, eta)
x1 = x(1);
x2 = x(2);
s = (2*pi*pi)*sin(pi*x1)*sin(pi*x2);
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
f = flux(u, q, w, v, x, t, mu, eta);
fb = f(1)*n(1) + f(2)*n(2) + tau*(u(1)-0.0);
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
x1 = x(1);
x2 = x(2);
ub = sin(pi*x1)*sin(pi*x2);
end

function u0 = initu(x, mu, eta)
u0 = sym(0.0);
end

function fb = fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau)
x1 = x(1);
x2 = x(2);
ub = sin(pi*x1)*sin(pi*x2);
fb1 = tau*(ub - uhat);
fb = [fb1 uhat];
end

function fi = fint(u, q, w, v, x, t, mu, eta, uhat, n, tau)
f = flux(u, q, w, v, x, t, mu, eta);
fi = f(1)*n(1) + f(2)*n(2) + tau*(u(1)-uhat(1));
end


