function flux(u, q, w, v, x, t, mu, eta)
    f = mu[1]*q + u[1]*mu[2:3] + u[1]*u[1]*mu[4:5];
    return f;
end
function source(u, q, w, v, x, t, mu, eta)
    s = 0.0;
    return s;
end
function ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    ub = [1-2*x[1] u];
    return ub;
end
function fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    f = flux(u, q, w, v, x, t, mu, eta);
    tm = f[1]*n[1] + f[2]*n[2] + tau[1]*(u[1]-uhat[1]);
    fb = [tm tm];
    return fb;
end
function initu(x, mu, eta)
    u0 = 1-2*x[1];
    return u0;
end
