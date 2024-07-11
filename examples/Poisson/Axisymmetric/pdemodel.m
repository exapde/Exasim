function pde = pdemodel
pde.flux = @flux;
pde.source = @source;
pde.fbou = @fbou;
pde.fbouhdg = @fbouhdg;
pde.ubou = @ubou;
pde.initu = @initu;
end

function f = flux(u, q, w, v, x, t, mu, eta)
z = x(1);
r = x(2);
f = mu*r*q;
end

function s = source(u, q, w, v, x, t, mu, eta)
z = x(1);
r = x(2);
s = sin(r)/exp(z);
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
f = flux(u, q, w, v, x, t, mu, eta);
tm = f(1)*n(1) + f(2)*n(2) + tau*(u(1)-uhat(1));
fb = [0.0 tm tm tm];
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
z = x(1);
r = x(2);
uexact = exp(-z)*cos(r);
ub = [u uexact uexact uexact]; 
end

function u0 = initu(x, mu, eta)
u0 = sym(0.0);
end

function fb = fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau)
f = flux(u, q, w, v, x, t, mu, eta);
fb1 = f(1)*n(1) + f(2)*n(2) + tau*(u(1)-uhat(1));
z = x(1);
r = x(2);
uexact = exp(-z)*cos(r);
fb2 = tau*(uexact - uhat);
fb = [fb1 fb2 fb2 fb2];
end


