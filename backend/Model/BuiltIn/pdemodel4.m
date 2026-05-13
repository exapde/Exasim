function pde = pdemodel4
    pde.mass = @mass;
    pde.flux = @flux;
    pde.source = @source;
    pde.fbouhdg = @fbouhdg;
    pde.fbou = @fbou;
    pde.ubou = @ubou;
    pde.initu = @initu;
    % pde.avfield = @avfield;
    pde.sourcew = @sourcew;
    pde.initw = @initw;
    pde.eos = @eos;
    pde.visscalars = @visscalars;
    pde.visvectors = @visvectors;
    pde.vistensors = @vistensors;
end

function m = mass(u, q, w, v, x, t, mu, eta)
    ns = 5;
    ndim = numel(x);
    m = sym(ones(ns + ndim + 1, 1));
end

function f = flux(u, q, w, v, x, t, mu, eta)
    f = fluxaxial2d(u, q, w, v, x, t, mu, eta);
end

function s = source(u, q, w, v, x, t, mu, eta)        
    s = sourceaxial2d(u, q, w, v, x, t, mu, eta);    
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    ns = 5;
    ndim = numel(x);
    ub = sym(ones(ns + ndim + 1, 1));
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    ns = 5;
    ndim = numel(x);
    fb = sym(ones(ns + ndim + 1, 1));
end

function u0 = initu(x, mu, eta)
    ns = 5;
    ndim = numel(x);
    u0 = sym(ones(ns + ndim + 1, 1));    
end

function w0 = initw(x, mu, eta)
    w0 = sym(ones(1,1));
end

function f = eos(u, q, w, v, x, t, mu, eta)
    f = eosnd(u, q, w, v, x, t, mu, eta);    
end

function f = sourcew(u, q, w, v, x, t, mu, eta)
    f = eosnd(u, q, w, v, x, t, mu, eta);        
end

function fb = fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau)
  fb = fbouhdgaxialnd(u, q, w, v, x, t, mu, eta, uhat, n, tau);
end

function s = visscalars(u, q, w, v, x, t, mu, eta)

[~, ~, ~, ~, rho_i_dim, T_dim, p_dim] = fluxaxial2d(u, q, w, v, x, t, mu, eta);

s = rho_i_dim; % physical (N, O, NO, N2, O2) densities
s(6) = sum(rho_i_dim); % physical mixture density
s(7) = T_dim; % phyiscal temperature
s(8) = p_dim; % physical pressure

end

function s = visvectors(u, q, w, v, x, t, mu, eta)

nd = length(x);

[~, ~, ~, ~, rho_i_dim, ~, ~, ~, heatflux_dim, chemflux_dim] = fluxaxial2d(u, q, w, v, x, t, mu, eta);

rho_scale   = mu(1);
u_scale     = mu(2);

rho = sum(rho_i_dim); % physical mixture density
s(:,1) = u(6:5+nd)*rho_scale*u_scale/rho; % physical velocity field
s(:,2) = heatflux_dim; % physical heat flux
s(:,3) = chemflux_dim; % physical chemical flux
s = s(:);

end

function s = vistensors(u, q, w, v, x, t, mu, eta)

[~, ~, ~, ~, ~, ~, ~, stress_dim, ~, ~] = fluxaxial2d(u, q, w, v, x, t, mu, eta);
% physical viscous stresses
s = stress_dim(:);

end

