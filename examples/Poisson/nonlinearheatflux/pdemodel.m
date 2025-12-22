
function pde = pdemodel
pde.flux = @flux;
pde.source = @source;
pde.sourcew = @sourcew;
pde.fbou = @fbou;
pde.ubou = @ubou;
pde.fbouhdg = @fbouhdg;
pde.initu = @initu;
pde.initw = @initw;
pde.eos = @eos;
end

function f = flux(u, q, w, v, x, t, mu, eta)
kappa = mu(1);
sw = sourcew(u, q, w, v, x, t, mu, eta);
dsdu  = jacobian(sw,u);
dsdw  = jacobian(sw,w);
dwdu = -dsdu/dsdw;
f = kappa*(dwdu*q);
end

function s = source(u, q, w, v, x, t, mu, eta)
x1 = x(1);
x2 = x(2);
s = x1*sin(x2);
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
f = flux(u, q, w, v, x, t, mu, eta);
fb = f(1)*n(1) + f(2)*n(2) + tau*(u(1)-0.0);
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
ub = sym(0.0); 
end

function u0 = initu(x, mu, eta)
u0 = sym(1.0);
end

function w0 = initw(x, mu, eta)
w0 = sym(0.693147180559945);
end

function s = eos(u, q, w, v, x, t, mu, eta)
s = exp(w) - sym(1.0) - u*u; 
end

function fb = fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau)
fb1 = tau*(mu(2) - uhat);
f = flux(uhat, q, w, v, x, t, mu, eta);
fb2 = f(1)*n(1) + f(2)*n(2) + tau*(u(1)-uhat(1));
fb = [fb1 fb2];
end

function s = sourcew(u, q, w, v, x, t, mu, eta)
s = exp(w(1)) - sym(1.0) - u*u; 
end


