function pde = pdemodel
pde.flux = @flux;
pde.source = @source;
pde.fbou = @fbou;
pde.fbouhdg = @fbouhdg;
pde.ubou = @ubou;
pde.initu = @initu;
pde.physparam = @param;
end

function f = flux(u, q, w, v, x, t, mu, eta)
r = x(1);
f = r*q;     % For axisymmetric
end

function s = source(u, q, w, v, x, t, mu, eta)    
    r_tilde = x(1);
    z_tilde = x(2);
    
    % Physics parameters
    l_ref = mu(2);
    E_ref = mu(4);
    e_eps0 = mu(5);
    N0 = mu(7);
    z0 = mu(8);
    sigma0 = mu(9);
    
    N0_tilde = N0*(l_ref^3);
    z0_tilde = z0/l_ref;
    sigma0_tilde = sigma0/l_ref;
    
    s = r_tilde.*(e_eps0/(E_ref*l_ref^2)).*N0_tilde.*exp(-((z_tilde - z0_tilde).^2 + r_tilde.^2)/(sigma0_tilde^2));
    
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
f = flux(u, q, w, v, x, t, mu, eta);
tm = f(1)*n(1) + f(2)*n(2) + tau*(u(1)-uhat(1));
fb = [tm sym(0) tm sym(0)]; % bottom, right, top, left
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)

l_ref = mu(2);
E_ref = mu(4);
e_eps0 = mu(5);
phi0 = mu(6);
N0 = mu(7);
z0 = mu(8);
sigma0 = mu(9);
phi0_tilde = phi0/(E_ref*l_ref);

ub = [sym(0) u phi0_tilde u]; % bottom, right, top, left
end

function u0 = initu(x, mu, eta)
u0 = sym(0.0);
end

function fb = fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau)
%fb = tau*(0.0 - uhat);

l_ref = mu(2);
E_ref = mu(4);
e_eps0 = mu(5);
phi0 = mu(6);
N0 = mu(7);
z0 = mu(8);
sigma0 = mu(9);
phi0_tilde = phi0/(E_ref*l_ref);

f = flux(u, q, w, v, x, t, mu, eta);
fb2 = f(1)*n(1) + f(2)*n(2) + tau*(u(1)-uhat(1));
fb1 = tau*(0 - uhat);
fb3 = tau*(phi0_tilde - uhat);
fb = [fb1 fb2 fb3 fb2];
end

