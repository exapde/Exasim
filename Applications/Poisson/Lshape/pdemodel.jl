
function flux(u, q, w, v, x, t, mu, eta)
    f = mu[1]*q;
    return f;
end
function source(u, q, w, v, x, t, mu, eta)
    s = 1.0;
    return s;
end
function ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    ub = 0.0;
    return ub;
end
function fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    f = flux(u, q, w, v, x, t, mu, eta);
    fb = f[1]*n[1] + f[2]*n[2] + tau[1]*(u[1]-0.0);
    return fb;
end
function initu(x, mu, eta)
    u0 = 0.0;
    return u0;
end
