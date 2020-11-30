
function pde = pdemodel2
pde.flux = @flux;
pde.source = @source;
pde.fbou = @fbou;
pde.ubou = @ubou;
pde.finterface = @finterface;
pde.uinterface = @uinterface;
pde.initu = @initu;
pde.initv = @initv;
end

function f = flux(u, q, w, v, x, t, mu, eta)
f = (mu+v(1))*q;
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
ub = sym(0.0); 
end

function fb = finterface(u, q, w, v, x, t, mu, eta, uhat, vhat, n, tau)
f = flux(u, q, w, v, x, t, mu, eta);
fb = f(1)*n(1) + f(2)*n(2) + tau*(u(1)-uhat(1));
end

function ub = uinterface(u, q, w, v, x, t, mu, eta, uhat, vhat, n, tau)
ub = vhat; 
end

function u0 = initu(x, mu, eta)
u0 = sym(0.0);
end

function v0 = initv(x, mu, eta)
v0 = sym(0.0);
end
