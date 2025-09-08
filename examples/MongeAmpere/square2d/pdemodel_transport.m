function pde = pdemodel
pde.mass = @mass;
pde.flux = @flux;
pde.source = @source;
pde.fbou = @fbou;
pde.fbouhdg = @fbouhdg;
pde.ubou = @ubou;
pde.initu = @initu;
end

function m = mass(u, q, w, v, x, t, mu, eta)
m = sym(1.0);
end

function f = flux(u, q, w, v, x, t, mu, eta)
v1 = v(1); % x-component of the velocity field
v2 = v(2); % y-component of the velocity field

theta = mu(2);
a1 = mu(3);
a2 = mu(4);
a = mu(5);
x1 = x(1);
x2 = x(2);

r = sqrt(x1^2 + x2^2);
rho = 1 + a1*sech(a2*(r^2 - a^2));
F = rho/theta;

w1 = v1/(t + (1-t)*F);
w2 = v2/(t + (1-t)*F);

f = [w1*u+mu(1)*(rho-1)*q(1) w2*u+mu(1)*(rho-1)*q(2)];

% F = v(1);
% v1 = v(4);
% v2 = v(5);
% w1 = v1/(t + (1-t)*F);
% w2 = v2/(t + (1-t)*F);
% f = [w1*u+mu(1)*(rho-1)*q(1) w2*u+mu(1)*(rho-1)*q(2)];

end

function s = source(u, q, w, v, x, t, mu, eta)

v1 = v(1);
v2 = v(2);

theta = mu(2);
a1 = mu(3);
a2 = mu(4);
a = mu(5);
x1 = x(1);
x2 = x(2);

r = sqrt(x1^2 + x2^2);
rho = 1 + a1*sech(a2*(r^2 - a^2));
rho1 = -(2*a1*a2*x1*sinh(a2*(- a^2 + x1^2 + x2^2)))/cosh(a2*(- a^2 + x1^2 + x2^2))^2;
rho2 = -(2*a1*a2*x2*sinh(a2*(- a^2 + x1^2 + x2^2)))/cosh(a2*(- a^2 + x1^2 + x2^2))^2;

F = rho/theta;
w1 = v1/(t + (1-t)*F);
w2 = v2/(t + (1-t)*F);

d = (w1*rho1 + w2*rho2)/theta;
s = u*(F - 1 - (1 - t)*d)/(t + (1-t)*F);

% F = v(1);
% F1 = v(2);
% F2 = v(3);
% v1 = v(4);
% v2 = v(5);
% w1 = v1/(t + (1-t)*F);
% w2 = v2/(t + (1-t)*F);
% d = (w1*F1 + w2*F2);
% s = u*(F - 1 - (1 - t)*d)/(t + (1-t)*F);

end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
f = flux(u, q, w, v, x, t, mu, eta);
fb = f(1)*n(1) + f(2)*n(2) + tau*(u(1)-uhat(1));
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
ub = sym(0.0); 
end

function u0 = initu(x, mu, eta)
u0 = x(1);
end

function fb = fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau)
fb = u(1) - uhat(1);
end


