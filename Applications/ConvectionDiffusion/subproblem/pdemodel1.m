
function pde = pdemodel1
pde.flux = @flux;
pde.source = @source;
pde.fbou = @fbou;
pde.ubou = @ubou;
pde.initu = @initu;
pde.initw = @initw;
end

function f = flux(u, q, w, v, x, t, mu, eta)
u1 = u(1);
q1x = q(1);
q1y = q(2);
q2x = w(2);
q2y = w(3);

fi = [u1*q2x u1*q2y];
fv = [mu(1)*q1x mu(1)*q1y];
f = fi + fv;
end

function s = source(u, q, w, v, x, t, mu, eta)
s = sym(1.0);
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

function w0 = initw(x, mu, eta)
w0 = [sym(0.0); sym(0.0); sym(0.0)];
end


