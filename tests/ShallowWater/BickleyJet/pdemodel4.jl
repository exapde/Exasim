function mass(u, q, w, v, x, t, mu, eta)
    m = [1.0, 1.0];
    return m;
end
function flux(u, q, w, v, x, t, mu, eta)
    r = w[1];  # pressure
    uv = u[1];
    vv = u[2];    
    p = r;    
    f = [uv*uv+p, vv*uv, uv*vv, vv*vv+p];
    f = reshape(f,(2,2));    
    return f;
end
function source(u, q, w, v, x, t, mu, eta)
    s = [0.0, 0.0];
    return s;
end
function sourcew(u, q, w, v, x, t, mu, eta)
    s = mu[1]*(q[1]+q[4]);
    return s;
end 
function ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    ub = [0.0, 0.0];
    return ub;
end
function fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    f = flux(u, q, w, v, x, t, mu, eta);
    fb = f[:,1]*n[1] + f[:,2]*n[2] + tau[1]*(u-uhat);
    return fb;
end
function initu(x, mu, eta)
    epsil = 0.1; # perturbation magnitude
    l = 0.5;     # Gaussian width
    k = 0.5;     # Sinusoidal wavenumber
    
    x1 = x[1];
    x2 = x[2];
    
    # The Bickley jet
    U = (1/cosh(x2))*(1/cosh(x2));

    # Slightly off-center vortical perturbations
    Psiprime = exp(-(x2 + l/10)*(x2 + l/10) / (2*(l*l))) * cos(k * x1) * cos(k * x2);

    # Vortical velocity fields (ũ, ṽ) = (-∂_y, +∂_x) ψ̃
    uprime =  Psiprime * (k * tan(k * x2) + x2 /(l*l)); 
    vprime = -Psiprime * k * tan(k * x1); 

    #u01 = 1.0; # h
    u01 = U + epsil * uprime; # h*u
    u02 = epsil * vprime;  # h*v
   
    u0 = [u01, u02];
    return u0;
end
function initw(x, mu, eta)
    w0 = [0.0];
    return w0;
end
