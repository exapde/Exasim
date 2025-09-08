function pde = pdemodel
pde.flux = @flux;
pde.source = @source;
%pde.sourcew = @sourcew;
pde.fbou = @fbou;
pde.fbouhdg = @fbouhdg;
pde.ubou = @ubou;
pde.initu = @initu;
%pde.initw = @initw;
end

function f = flux(u, q, w, v, x, t, mu, eta)
Q = reshape(q, [2 2]);
f = mu(1)*(Q + Q.') + mu(2)*(Q(1,1) + Q(2,2))*eye(2,2);
end

function s = source(u, q, w, v, x, t, mu, eta)
s = 0*x;
nu = mu(1);
x1 = x(1);
x2 = x(2);
Re = 1/nu; 
lam = Re/2 - sqrt(Re^2+16*pi^2)/2;
s(1) = lam*exp(lam*x1)*(cos(2*pi*x2) - exp(lam*x1));
s(2) = -(lam^2*exp(lam*x1)*sin(2*pi*x2))/(2*pi);
end

function s = sourcew(u, q, w, v, x, t, mu, eta)
s = w(1) - mu(2)*(q(1)+q(4));
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

function w0 = initw(x, mu, eta)
w0 = sym(0.0);
end

function fb = fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau)
Re = 1/(mu(1)); 
lam = Re/2 - sqrt(Re^2/4+4*pi^2);
x1 = x(1);
x2 = x(2);
ue = 0*x;
ue(1) = 1-exp(lam*x1)*cos(2*pi*x2);
ue(2) = (lam/(2*pi))*exp(lam*x1)*sin(2*pi*x2);
fb = tau*(ue - uhat(:));
end

% function [u,v,p] = exacsol(mu, x, y) 
%     Re = 1/(mu(1)); 
%     lam = Re/2 - sqrt(Re^2+16*pi^2)/2;
%     u = 1-exp(lam*x)*cos(2*pi*y);
%     v = (lam/(2*pi))*exp(lam*x)*sin(2*pi*y);
%     p = 0.5 - 0.5*exp(2*lam*x);
%     ux = -lam*exp(lam*x)*cos(2*pi*y);
%     uy = 2*pi*exp(lam*x)*sin(2*pi*y);
%     vx = (lam^2*exp(lam*x)*sin(2*pi*y))/(2*pi);
%     vy = lam*exp(lam*x)*cos(2*pi*y);
%     px = -lam*exp(2*lam*x);
%     py = 0;
%     uxx = -lam^2*exp(lam*x)*cos(2*pi*y);
%     uyy = 4*pi^2*exp(lam*x)*cos(2*pi*y);
%     vxx = (lam^3*exp(lam*x)*sin(2*pi*y))/(2*pi);
%     vyy = -2*lam*pi*exp(lam*x)*sin(2*pi*y);
% 
%     sx = simplify(-mu(1)*(uxx + uyy) + px + u*ux + v*uy);
%     sy = simplify(-mu(1)*(vxx + vyy) + py + u*vx + v*vy);    
%     fx = simplify(-u*ux - v*uy);
%     fy = simplify(-u*vx - v*vy);    
% end
 
