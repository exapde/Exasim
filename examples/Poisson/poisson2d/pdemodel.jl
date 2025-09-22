
function flux(u, q, w, v, x, t, mu, eta)
    f = mu[1]*q;
    return f;
end
function source(u, q, w, v, x, t, mu, eta)
    x1 = x[1];
    x2 = x[2];
    s = (2*pi*pi)*sin(pi*x1)*sin(pi*x2);
    return s;
end
function visscalars(u, q, w, v, x, t, mu, eta)
    s = 0*q;
    s[1] = u[1];
    s[2] = q[1] + q[2];
    return s;
end
function visvectors(u, q, w, v, x, t, mu, eta)
    return q
end
function qoivolume(u, q, w, v, x, t, mu, eta)
    s = 0*q;
    x1 = x[1];
    x2 = x[2];
    uexact = sin(pi*x1)*sin(pi*x2);
    s[1] = (u[1] - uexact)^2;
    s[2] = u[1];
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
function qoiboundary(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    f = flux(u, q, w, v, x, t, mu, eta);
    fb = f[1]*n[1] + f[2]*n[2] + tau[1]*(u[1]-uhat[1]);
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
