function flux(u, q, w, v, x, t, mu, eta)
    f = mu[1]*q + u[1]*mu[2];
    return f;
end
function source(u, q, w, v, x, t, mu, eta)
    x1 = x[1];
    nu = mu[1];
    c = mu[2];
    c2 = c^2;
    s = c2*nu*(sin(c*x1) - cos(c*x1));
    return s;
end
function ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    ub = 0.0;
    return ub;
end
function fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    f = flux(u, q, w, v, x, t, mu, eta);
    fb = f[1]*n[1] + tau[1]*(u[1]-uhat[1]);
    return fb;
end
function initu(x, mu, eta)
    u0 = 0.0;
    return u0;
end
