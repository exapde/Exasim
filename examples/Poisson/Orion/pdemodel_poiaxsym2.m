function pde = pdemodel
% -(1/r) * d/dr (r du/dr) - d2u/dz2 = s
% -(1/r) du/dr - d2u/dr2 - d2u/dz2 = s
% - d2u/dr2 - d2u/dz2 = s + (1/r) du/dr

pde.flux = @flux;
pde.source = @source;
pde.fbou = @fbou;
pde.fbouhdg = @fbouhdg;
pde.ubou = @ubou;
pde.initu = @initu;
end

function f = flux(u, q, w, v, x, t, mu, eta)
f = mu(1)*q;
end

function s = source(u, q, w, v, x, t, mu, eta)
r = x(2);
s = mu(2) - mu(1)*q(2)/r;
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
fb = sym(0.0);
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
ub = sym(0.0);
end

function u0 = initu(x, mu, eta)
u0 = sym(0.0);
end

function fb = fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau)
f = flux(u, q, w, v, x, t, mu, eta);
fb1 = f(1)*n(1) + f(2)*n(2) + tau*(u(1)-uhat(1));
fb2 = tau*(0.0 - uhat);
fb3 = tau*(1.0 - uhat);
fb4 = fb1 + mu(4)*(mu(3) - uhat(1));
fb = [fb1 fb2 fb3 fb4];
end


