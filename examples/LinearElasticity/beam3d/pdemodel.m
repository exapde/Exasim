function pde = pdemodel
pde.flux = @flux;
pde.source = @source;
pde.fbou = @fbou;
pde.fbouhdg = @fbouhdg;
pde.ubou = @ubou;
pde.initu = @initu;
end

function f = flux(u, q, w, v, x, t, mu, eta)
Q = reshape(q, [3 3]);
f = mu(1)*(Q + Q.') + mu(2)*(Q(1,1) + Q(2,2) + Q(3,3))*eye(3,3);
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
f = flux(u, q, w, v, x, t, mu, eta);
fh = f(:,1)*n(1) + f(:,2)*n(2) + f(:,3)*n(3) + tau*(u(:)-uhat(:));
fb2 = tau*(0*x - uhat(:));
fi = mu(3:5);
fb3 = fh - fi(:);
fb = [fh fb2 fb3];
end

