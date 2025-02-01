
function pde = pdemodel
pde.flux = @flux;
pde.source = @source;
pde.fbou = @fbou;
pde.fbouhdg = @fbouhdg;
pde.ubou = @ubou;
pde.initu = @initu;
end

function f = flux(u, q, w, v, x, t, mu, eta)
f = mu*q;
end

function s = source(u, q, w, v, x, t, mu, eta)
x1 = x(1);

%s = 1*x1/x1; %examples 1 and 2
s = -4*(x1-0.5)^2+1; %example 3 and 4
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
f = flux(u, q, w, v, x, t, mu, eta);

fb(1) = sym(0.0);
fb(2) = mu*q;
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
ub(1) = sym(0.0);
ub(2) = u(1);
end

function u0 = initu(x, mu, eta)
x1 = x(1);
u0 = 0*(x1*x1);
end

function fb = fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau)
fb = tau*(x^4/3-2*x^3/3 - uhat);
end


