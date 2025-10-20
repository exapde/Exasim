function pde = pdemodel
pde.flux = @flux;
pde.source = @source;
pde.fbou = @fbou;
pde.fbouhdg = @fbouhdg;
pde.ubou = @ubou;
pde.initu = @initu;
end

function f = flux(u, q, w, v, x, t, mu, eta)

kappa = mu(1);
lambda = mu(2);

F11 = q(1);
F21 = q(2);
F31 = q(3);
F12 = q(4);
F22 = q(5);
F32 = q(6);
F13 = q(7);
F23 = q(8);
F33 = q(9);

F = -[F11 F12 F13; F21 F22 F23; F31 F32 F33];
JF = det(F);
CF = F.'*F;
W  = (kappa/2)*(trace(CF)-3) - kappa*log(JF) + (lambda/2)*(JF-1)^2;
f = reshape(jacobian(W,q(:)), [3 3]);

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
fh = f(:,1)*n(1) + f(:,2)*n(2) + f(:,3)*n(3) + tau*(u(:)-uhat(:));
fb2 = (x - uhat(:));
fi = mu(3:5);
fb3 = fh + fi(:);
fb = [fh fb2 fb3];
end

