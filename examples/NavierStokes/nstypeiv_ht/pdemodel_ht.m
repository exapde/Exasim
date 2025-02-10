
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
f = flux(u, q, w, v, x, t, mu, eta);
fh = f(1)*n(1) + f(2)*n(2) + tau*(u(1)-uhat(1));
fb1 = tau*(mu(2) - uhat(1));
fb2 = tau*(v(1) - uhat(1));
fb3 = fh;
fh_ext = -(v(2)*n(1) + v(3)*n(2)) + tau*(v(1)-v(4));
fb4 = fh + fh_ext;
% Direchlet, Direchlet, Neumann, Neumann
fb = [fb1; fb2; fb3; fb4];
end




