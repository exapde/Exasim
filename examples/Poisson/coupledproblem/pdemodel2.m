
function pde = pdemodel2
pde.flux = @flux;
pde.source = @source;
pde.fbou = @fbou;
pde.fint = @fint;
pde.fbouhdg = @fbouhdg;
pde.ubou = @ubou;
pde.initu = @initu;
end

function f = flux(u, q, w, v, x, t, mu, eta)
f = mu*q;
end

function s = source(u, q, w, v, x, t, mu, eta)
x1 = x(1);
x2 = x(2);
s = (2*pi*pi)*sin(pi*x1)*sin(pi*x2);
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
f = flux(u, q, w, v, x, t, mu, eta);
fb = f(1)*n(1) + f(2)*n(2) + tau*(u(1)-0.0);
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
ub = sym(0.0); 
end

function u0 = initu(x, mu, eta)
x1 = x(1);
x2 = x(2);
u0 = 0*x1*x2*(1-x1)*(1-x2);
%u0 = sin(pi*x1)*sin(pi*x2);
end

function fb = fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau)
x1 = x(1);
x2 = x(2);
ub = sin(pi*x1)*sin(pi*x2);
qx = -pi*cos(pi*x1)*sin(pi*x2);
qy = -pi*sin(pi*x1)*cos(pi*x2);
fb1 = (0 - uhat);
fb = [fb1 uhat-0*ub];
% f = flux(u, q, w, v, x, t, mu, eta);
% fi = f(1)*n(1) + f(2)*n(2) + tau*(u(1)-uhat(1));
% fb = [fb1 fi-0*(qx*n(1)+qy*n(2))];
end

function fi = fint(u, q, w, v, x, t, mu, eta, uhat, n, tau)
%fi = -uhat(1);
f = flux(u, q, w, v, x, t, mu, eta);
fi =  ( f(1)*n(1) + f(2)*n(2) + tau*(u(1)-uhat(1)) );
end


