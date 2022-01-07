function mass(u, q, w, v, x, t, mu, eta)
    m = 1.0;
    return m;
end
function flux(u, q, w, v, x, t, mu, eta)
    f = mu[1]*q;
    return f;
end
function source(u, q, w, v, x, t, mu, eta)    
    s = 0.0;
    return s;
end
function ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    ub = mu[2];
    return ub;
end
function fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    f = flux(u, q, w, v, x, t, mu, eta);
    fb = f[1]*n[1] + f[2]*n[2] + tau[1]*(u[1]-0.0);
    return fb;
end
function initu(x, mu, eta)
    x1 = x[1];
    x2 = x[2];
    u0 = sin(pi*x1)*sin(pi*x2);    
    return u0;
end
function initq(x, mu, eta)
    q0 = [0.0, 0.0];
end
function initw(x, mu, eta)
    w0 = 0.0;
    return w0;
end
