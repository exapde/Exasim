
function pde = pdemodel
pde.flux = @flux;
pde.source = @source;
pde.fbou = @fbou;
pde.ubou = @ubou;
pde.initu = @initu;
end

function f = flux(u, q, w, v, x, t, mu, eta)
f = mu(1)*q + u.^2 / 2;
end

function s = source(u, q, w, v, xdg, t, mu, eta)
x = xdg(1);
y = xdg(2);
k = mu(1);

wx = (1-x)/k;
wy = (1-y)/k;

u = x*y*tanh(wx)*tanh(wy);
u_x = (-x/k * (sech(wx))^2 + tanh(wx)) * y * tanh(wy);
u_xx = -2*(sech(wx))^2 / k^2 * (x * tanh(wx) + k) * y * tanh(wy);
u_y = (-y/k * (sech(wy))^2 + tanh(wy)) * x * tanh(wx);
u_yy = -2*(sech(wy))^2 / k^2 * (y * tanh(wy) + k) * x * tanh(wx);

s = -k * u_xx - k * u_yy + u * u_x + u * u_y;
% 
% s = x*y*((x*y*(tanh((1 - x) / k)^2 - 1)*tanh((1 - y) / k)) / k + y*tanh((1 - x) / k)*tanh((1 - y) / k))*tanh((1 - x) / k)*tanh((1 - y) / k) ...
%     + x*y*((x*y*(tanh((1 - y) / k)^2 - 1)*tanh((1 - x) / k)) / k + x*tanh((1 - x) / k)*tanh((1 - y) / k))*tanh((1 - x) / k)*tanh((1 - y) / k) ...
%     - k*(((2x*y*(tanh((1 - x) / k)^2 - 1)*tanh((1 - x) / k)*tanh((1 - y) / k)) / k ...
%     + y*(tanh((1 - x) / k)^2 - 1)*tanh((1 - y) / k)) / k ...
%     + (y*(tanh((1 - x) / k)^2 - 1)*tanh((1 - y) / k)) / k) ...
%     - k*(((2x*y*(tanh((1 - y) / k)^2 - 1)*tanh((1 - x) / k)*tanh((1 - y) / k)) / k ...
%     + x*(tanh((1 - y) / k)^2 - 1)*tanh((1 - x) / k)) / k + (x*(tanh((1 - y) / k)^2 - 1)*tanh((1 - x) / k)) / k)
% s = sym(0.0); 
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
f = flux(u, q, w, v, x, t, mu, eta);
tm = f(1)*n(1) + f(2)*n(2) + tau*(u(1)-uhat(1));
fb = [tm, tm];
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
ub = [0*u, 0*u]; 
end

function u0 = initu(x, mu, eta)
u0 = sym(0.0);
end
