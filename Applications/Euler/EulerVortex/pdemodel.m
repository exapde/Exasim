% Poisson equation with homogenous Dirichlet condition on a unit square
% q + \nabla u = 0 in \Omega
% \nabla dot flux = source in \Omega
% u = 0 on \partial Omega
% flux = 2*param*q
% source = 2*pi*pi*sin(pi*x1)*sin(pi*x2);

function pde = pdemodel
pde.mass = @mass;
pde.flux = @flux;
pde.source = @source;
pde.fbou = @fbou;
pde.ubou = @ubou;
pde.initu = @initu;
end

function m = mass(u, q, w, v, x, t, mu, eta)
m = sym([1.0; 1.0; 1.0; 1.0]); 
end

function f = flux(u, q, w, v, x, t, mu, eta)
    gam = mu(1);
    gam1 = gam - 1.0;
    r = u(1);
    ru = u(2);
    rv = u(3);
    rE = u(4);
    r1 = 1/r;
    uv = ru*r1;
    vv = rv*r1;
    E = rE*r1;
    p = gam1*(rE-r*0.5*(uv*uv+vv*vv));
    h = E+p*r1;
    f = [ru, ru*uv+p, rv*uv, ru*h, rv, ru*vv, rv*vv+p, rv*h];
    f = reshape(f,[4,2]);    
end

function s = source(u, q, w, v, x, t, mu, eta)
s = [sym(0.0); sym(0.0); sym(0.0); sym(0.0)];
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
f = flux(u, q, w, v, x, t, mu, eta);
fb = f(:,1)*n(1) + f(:,2)*n(2) + tau*(u-uhat);
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
ub = sym([0.0; 0.0; 0.0; 0.0]); 
end

function u0 = initu(x, mu, eta)

    t = 0;
    phi = 5;
    x1 = x(1);
    x2 = x(2);
    gam = mu(1);
    M_ref = mu(2);
    
    r = sqrt((x1-t).^2 + x2.^2);
    u01 = (1 - ((gam-1)/(16*pi^2)) * phi^2 * exp(2*(1-r.^2)) ).^(1/(gam-1));
    u02 = u01 .* (1 - M_ref^(-1) * phi * exp(1-r.^2) .* x2/(2*pi));
    u03 = u01 .* M_ref^(-1) * phi .* exp(1-r.^2) .* x1/(2*pi);
    p = u01.^gam/(gam*M_ref^2);
    u04 = p/(gam-1) + 0.5 * (u02.*u02 + u03.*u03) ./ u01;

    u0 = [u01; u02; u03; u04];
end

