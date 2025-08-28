function pde = pdemodel
pde.flux = @flux;
pde.source = @source;
pde.fbou = @fbou;
pde.fbouhdg = @fbouhdg;
pde.ubou = @ubou;
pde.initu = @initu;
end

function f = flux(u, q, w, v, x, t, mu, eta)
F11 = q(1);
F21 = q(2);
F12 = q(3);
F22 = q(4);
F=[F11 F12; F21 F22];
detF = det(F);
lambda = mu(2);
f = [F11*mu(1) - F22*lambda*(F12*F21 - F11*F22 + 1) - (F22*mu(1))/(detF) ...
     F12*mu(1) + F21*lambda*(F12*F21 - F11*F22 + 1) + (F21*mu(1))/(detF) ; ...
     F21*mu(1) + F12*lambda*(F12*F21 - F11*F22 + 1) + (F12*mu(1))/(detF) ...
     F22*mu(1) - F11*lambda*(F12*F21 - F11*F22 + 1) - (F11*mu(1))/(detF)];
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
u0 = x;
end

function fb = fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau)
f = flux(u, q, w, v, x, t, mu, eta);
fh = f(:,1)*n(1) + f(:,2)*n(2) + tau*(u(:)-uhat(:));
fb2 = tau*(x - uhat(:));
fi = mu(3:4);
fb3 = fh + fi(:);
fb = [fh fb2 fb3];
end

