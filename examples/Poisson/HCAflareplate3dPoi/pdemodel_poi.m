function pde = pdemodel
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
s = sym(0.0);
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
fb1 = (0.0 - uhat);
fb2 = (1.0 - uhat);
f = flux(uhat, q, w, v, x, t, mu, eta);
fb3 = f(1)*n(1) + f(2)*n(2) + f(3)*n(3) +  tau*(u(1)-uhat(1));
fb = [fb1 fb2 fb3];
end


