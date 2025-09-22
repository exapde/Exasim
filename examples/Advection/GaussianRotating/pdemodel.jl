function mass(u, q, w, v, x, t, mu, eta)
    m = 1.0;
    return m;
end
function flux(u, q, w, v, x, t, mu, eta)
    f = [mu[1]*x[2]*u[1] -mu[2]*x[1]*u[1]];
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
    u0 = exp(-(x[1]-0.25)*(x[1]-0.25)/0.01-(x[2]-0.0)*(x[2]-0.0)/0.01);
    return u0;
end

