
function pde = pdemodel
pde.flux = @flux;
pde.source = @source;
pde.fbou = @fbou;
pde.fbouhdg = @fbouhdg;
pde.ubou = @ubou;
pde.initu = @initu;
pde.visscalars = @visscalars;
pde.visvectors = @visvectors;
pde.qoivolume = @qoivolume;
pde.qoiboundary = @qoiboundary;
end

function f = flux(u, q, w, v, x, t, mu, eta)
f = mu*q;
end

function s = source(u, q, w, v, x, t, mu, eta)
x1 = x(1);
x2 = x(2);
s = (2*pi*pi)*sin(pi*x1)*sin(pi*x2);
end

function s = visscalars(u, q, w, v, x, t, mu, eta)
s(1) = u(1);
s(2) = q(1) + q(2);
end

function s = visvectors(u, q, w, v, x, t, mu, eta)
s = q;
end

function s = qoivolume(u, q, w, v, x, t, mu, eta)
x1 = x(1);
x2 = x(2);
uexact = sin(pi*x1)*sin(pi*x2);
s(1) = (u(1) - uexact)^2;
s(2) = u(1);
end

function sb = qoiboundary(u, q, w, v, x, t, mu, eta, uhat, n, tau)
f = flux(u, q, w, v, x, t, mu, eta);
sb = f(1)*n(1) + f(2)*n(2) + tau*(u(1)-uhat(1));
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
u0 = 0*(x1*x1 + x2*x2);
end

function fb = fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau)
fb = tau*(0.0 - uhat);
end


