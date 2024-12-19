function pde = pdemodel
pde.flux = @flux;
pde.source = @source;
pde.fbouhdg = @fbouhdg;
pde.initu = @initu;
end

function f = flux(u, q, w, v, x, t, mu, eta)
% - nabla dot (kappa nabla w(u)) = s
% w(u) = log(1 + u^2)
% dwdu = (2*u)/(u^2 + 1)
% nabla w(u) = dwdu * nabla u
kappa = mu(1);
w = log(u^2 + 1);
dwdu = (2*u)/(u^2 + 1);
f = kappa*(dwdu*q);
end

function s = source(u, q, w, v, x, t, mu, eta)
x1 = x(1);
x2 = x(2);
s = x1*sin(x2);
end

function u0 = initu(x, mu, eta)
u0 = sym(1.0);
end

function fb = fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau)
fb1 = tau*(mu(2) - uhat);
f = flux(uhat, q, w, v, x, t, mu, eta);
fb2 = f(1)*n(1) + f(2)*n(2) + tau*(u(1)-uhat(1));
fb = [fb1 fb2];
end


