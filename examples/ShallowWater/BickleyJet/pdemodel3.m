function pde = pdemodel3
pde.mass = @mass;
pde.flux = @flux;
pde.source = @source;
pde.sourcew = @sourcew;
pde.fbou = @fbou;
pde.ubou = @ubou;
pde.initu = @initu;
pde.initw = @initw;
end

function m = mass(u, q, w, v, x, t, mu, eta)
m = sym([1.0; 1.0]); 
end

function f = flux(u, q, w, v, x, t, mu, eta)
    g = mu(1);    
    r = w(1);  % total height: h
    ru = u(1); % h*u 
    rv = u(2); % h*v   
    r1 = 1/r;
    uv = ru*r1; %u
    vv = rv*r1; %v    
    p = (0.5*g)*(r*r); % 0.5*g*h^2    
    f = [ru*uv+p, rv*uv, ru*vv, rv*vv+p];
    f = reshape(f,[2,2]);    
end

function s = source(u, q, w, v, x, t, mu, eta)
s = [sym(0.0); sym(0.0)];
end

function s = sourcew(u, q, w, v, x, t, mu, eta)
s = q(1)+q(4);
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
f = flux(u, q, w, v, x, t, mu, eta);
fb1 = f(:,1)*n(1) + f(:,2)*n(2) + tau*(u-uhat);
fw2 = f(:,1)*n(1) + f(:,2)*n(2) + tau*(u-uhat);
g = mu(1);   
r = w(1);
p = (0.5*g)*(r*r);
fw2(1) = p*n(1);
fw2(2) = p*n(2);
fb = [fb1 fw2];
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
ub = sym([0.0 0.0; 0.0 0.0]); 
end

function u0 = initu(x, mu, eta)

    epsil = 0.1; % perturbation magnitude
    l = 0.5;     % Gaussian width
    k = 0.5;     % Sinusoidal wavenumber
    
    x1 = x(1);
    x2 = x(2);
    
    % The Bickley jet
    U = sech(x2)*sech(x2);

    % Slightly off-center vortical perturbations
    Psiprime = exp(-(x2 + l/10)*(x2 + l/10) / (2*(l*l))) * cos(k * x1) * cos(k * x2);

    % Vortical velocity fields (ũ, ṽ) = (-∂_y, +∂_x) ψ̃
    uprime =  Psiprime * (k * tan(k * x2) + x2 /(l*l)); 
    vprime = -Psiprime * k * tan(k * x1); 

    %u01 = 1.0; % h
    u01 = U + epsil * uprime; % h*u
    u02 = epsil * vprime;  % h*v
    
    u0 = [u01; u02];
end

function w0 = initw(x, mu, eta)
w0 = sym(1.0);
end

