
function pde = pdemodel
pde.flux = @flux;
pde.source = @source;
pde.fbou = @fbou;
pde.fbouhdg = @fbouhdg;
pde.ubou = @ubou;
pde.initu = @initu;
end

function f = flux(u, q, w, v, x, t, mu, eta)
f = mu(1)*q + u(1)*mu(2:3) + u(1)*u(1)*mu(4:5);
end

function s = source(u, q, w, v, x, t, mu, eta)
s = sym(0.0); 
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
f = flux(u, q, w, v, x, t, mu, eta);
tm = f(1)*n(1) + f(2)*n(2) + tau*(u(1)-uhat(1));
fb = [tm, tm];
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
ub = [1-2*x(1), u]; 
end

function u0 = initu(x, mu, eta)
u0 = 1 - 2*x(1);
end

function fb = fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau)
fb1 = tau*(1-2*x(1)-uhat(1));
fb2 = tau*(u(1)-uhat(1));
fb = [fb1 fb2];
end



