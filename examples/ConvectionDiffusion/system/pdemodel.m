
function pde = pdemodel
pde.flux = @flux;
pde.source = @source;
pde.fbou = @fbou;
pde.fbouhdg = @fbouhdg;
pde.ubou = @ubou;
pde.initu = @initu;
end

function f = flux(u, q, w, v, x, t, mu, eta)
u1 = u(1);
u2 = u(2);
q1x = q(1);
q2x = q(2);
q1y = q(3);
q2y = q(4);
fi = [u1*q2x sym(0) u1*q2y sym(0)];
fv = [mu(1)*q1x q2x mu(1)*q1y q2y];
f = fi + fv;
f = reshape(f,[2,2]);    
end

function s = source(u, q, w, v, x, t, mu, eta)
u1 = u(1);
s = [sym(1.0); u1]; 
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
f = flux(u, q, w, v, x, t, mu, eta);
fb = f(:,1)*n(1) + f(:,2)*n(2) + tau*(u-uhat);
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
ub = [sym(0.0); sym(0.0)]; 
end

function u0 = initu(x, mu, eta)
u0 = [sym(0.0); sym(0.0)];
end

function fb = fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau)
fb0 = tau*(sym(0.0)-uhat(1));
fb = [fb0; fb0];
end



