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
f = mu(1) * (Q + Q.') + mu(2) * (Q(1,1) + Q(2,2)) * eye(2,2);
end

function s = source(u, q, w, v, x, t, mu, eta)
s = 0 * x;
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
fb = 0 * x;
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
ub = 0 * x;
end

function u0 = initu(x, mu, eta)
u0 = 0 * x;
end

function fb = fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau)
% f = flux(uhat, q, w, v, x, t, mu, eta);
% fh = f(:,1) * n(1) + f(:,2) * n(2) + tau * (u(:) - uhat(:));
ub = 0*uhat;
amp = mu(3);
ub(1) = 0.0;
ub(2) = amp * sin(2 * pi * x(1));

fb1 = tau * (0*ub(:) - uhat(:));
fb2 = tau * (ub(:) - uhat(:));
fb = [fb1 fb2];
end
