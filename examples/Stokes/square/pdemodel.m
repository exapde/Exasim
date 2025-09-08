function pde = pdemodel
pde.flux = @flux;
pde.source = @source;
pde.fbou = @fbou;
pde.fbouhdg = @fbouhdg;
pde.ubou = @ubou;
pde.initu = @initu;
end

function f = flux(u, q, w, v, x, t, mu, eta)
Q = reshape(q, [2 2]);

theta = mu(3);
a1 = mu(4);
a2 = mu(5);
a = mu(6);
x1 = x(1);
x2 = x(2);
r = sqrt(x1^2 + x2^2);
rho = 1 + a1*sech(a2*(r^2 - a^2));
F = rho/theta;
s = 1 - F;

f = mu(1)*(Q + 0*Q.') + mu(2)*(Q(1,1) + Q(2,2) - s)*eye(2,2);
end

function s = source(u, q, w, v, x, t, mu, eta)
s = 0*x;
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
fb = 0*x;
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
ub = 0*x; 
end

function u0 = initu(x, mu, eta)
u0 = 0*x;
end

function fb = fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau)
ue = 0*x;
fb = tau*(ue - uhat(:));
end

