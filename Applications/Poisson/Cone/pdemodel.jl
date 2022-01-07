function flux(u, q, w, v, x, t, mu, eta)
    f = mu[1]*q;
    return f;
end
function source(u, q, w, v, x, t, mu, eta)
    s = 0.0;
    return s;
end
function ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    ub = [mu[2] u];
    return ub;
end
function fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    f = flux(u, q, w, v, x, t, mu, eta);
    tm = f[1]*n[1] + f[2]*n[2] + f[3]*n[3] + tau[1]*(u[1]-uhat[1]);
    fb = [tm mu[3]*n[1]+mu[4]*n[2]+mu[5]*n[3]];
    return fb;
end
function initu(x, mu, eta)
    u0 = 1.0;
    return u0;
end
