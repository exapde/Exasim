
function pde = pdemodel
pde.flux = @flux;
pde.source = @source;
pde.fbou = @fbou;
pde.fbouhdg = @fbouhdg;
pde.ubou = @ubou;
pde.initu = @initu;
end

function f = flux(u, q, w, v, x, t, mu, eta)
f = mu*v(1)*q;
end

function s = source(u, q, w, v, x, t, mu, eta)
s = v(2) - u;
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
f = flux(u, q, w, v, x, t, mu, eta);
fb = f(1)*n(1) + f(2)*n(2) + tau*(u(1)-0.0);
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
ub = sym(0.0); 
end

function u0 = initu(x, mu, eta)
u0 = sym(0.0);
end

function fb = fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau)
f = flux(u, q, w, v, x, t, mu, eta);
fb1 = f(1)*n(1) + f(2)*n(2) + tau*(u-uhat);
fb2 = tau*(0-uhat); 
fb3 = tau*(v(2)-uhat); 
fb4 = f(1)*n(1) + f(2)*n(2) + tau*(u-uhat) + tau*(v(2)-uhat);
fb = [fb1 fb2 fb3 fb4];
end



