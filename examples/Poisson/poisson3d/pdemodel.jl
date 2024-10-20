function flux(u, q, w, v, x, t, mu, eta)
    f = mu[1]*q;
    return f;
end
function source(u, q, w, v, x, t, mu, eta)
    s = (3*pi*pi)*sin(pi*x[1])*sin(pi*x[2])*sin(pi*x[3]);
    return s;
end
function ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    ub = mu[2];
    return ub;
end
function fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    f = flux(u, q, w, v, x, t, mu, eta);
    fb = f[1]*n[1] + f[2]*n[2] + f[3]*n[3] + tau[1]*(u[1]-0.0);
    return fb;
end
function fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    fb = tau[1]*(0.0 - uhat[1]);
    return fb;
end
function initu(x, mu, eta)
    u0 = 0.0;
    return u0;
end
