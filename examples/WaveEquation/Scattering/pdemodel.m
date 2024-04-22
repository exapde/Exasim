function pde = pdemodel
pde.mass = @mass;
pde.flux = @flux;
pde.source = @source;
pde.fbou = @fbou;
pde.ubou = @ubou;
pde.initu = @initu;
pde.initq = @initq;
pde.initw = @initw;
end

function m = mass(u, q, w, v, x, t, mu, eta)
m = sym(1.0);
end

function f = flux(u, q, w, v, x, t, mu, eta)
f = mu(1)*q;
end

function s = source(u, q, w, v, x, t, mu, eta)
s = sym(0.0);
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
%ui = sin(mu(2)*x(1)+mu(3)*x(2)-mu(4)*sqrt(mu(2)*mu(2)+mu(3)*mu(3))*t);
uix = mu(2)*cos(mu(2)*x(1)+mu(3)*x(2)-mu(4)*sqrt(mu(2)*mu(2)+mu(3)*mu(3))*t);
uiy = mu(3)*cos(mu(2)*x(1)+mu(3)*x(2)-mu(4)*sqrt(mu(2)*mu(2)+mu(3)*mu(3))*t);
fb = [u uix*n(1)+uiy*n(2)]; % Neumann conditions
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
ub = [u u]; % Neumann conditions
end

function u0 = initu(x, mu, eta)
u0 = sym(0.0);
end

function q0 = initq(x, mu, eta)
q0 = [sym(0.0) sym(0.0)];
end

function w0 = initw(x, mu, eta)
w0 = sym(0.0);
end


