
function pde = pdemodel4
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
    %g = mu(1);        
    r = w(1);  % pressure
    uv = u(1);
    vv = u(2);    
    p = r;    
    f = [uv*uv+p, vv*uv, uv*vv, vv*vv+p];
    f = reshape(f,[2,2]);    
end

function s = source(u, q, w, v, x, t, mu, eta)
s = [sym(0.0); sym(0.0)];
end

function s = sourcew(u, q, w, v, x, t, mu, eta)
g = mu(1);    
s = g*(q(1)+q(4));
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
f = flux(u, q, w, v, x, t, mu, eta);
fb1 = f(:,1)*n(1) + f(:,2)*n(2) + tau*(u-uhat);
fw2 = f(:,1)*n(1) + f(:,2)*n(2) + tau*(u-uhat);
%g = mu(1);   
r = w(1);
p = r;
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

    %u01 = 0.0;
    u01 = U + epsil * uprime;
    u02 = epsil * vprime;
    
    u0 = [u01; u02];
end

function w0 = initw(x, mu, eta)
w0 = sym(0.0);
end


