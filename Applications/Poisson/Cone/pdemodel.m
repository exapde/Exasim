function pde = pdemodel
pde.flux = @flux;
pde.source = @source;
pde.fbou = @fbou;
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
f = flux(u, q, w, v, x, t, mu, eta);
tm = f(1)*n(1) + f(2)*n(2) + f(3)*n(3) + tau*(u(1)-uhat(1));
fb = [tm mu(3)*n(1)+mu(4)*n(2)+mu(5)*n(3)];
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
ub = [mu(2) u]; 
end

function u0 = initu(x, mu, eta)
u0 = sym(1.0);
end


