
function flux(u, q, w, v, x, t, mu, eta)
    z = x[1];
    r = x[2];
    f = mu[1]*r*q;
    return f;
end
function source(u, q, w, v, x, t, mu, eta)
    z = x[1];
    r = x[2];
    s = sin(r)/exp(z);
    return s;
end
function ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    z = x[1];
    r = x[2];
    uexact = exp(-z)*cos(r);
    ub = [u uexact uexact uexact]; 
    return ub;
end
function fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    f = flux(u, q, w, v, x, t, mu, eta);
    tm = f[1]*n[1] + f[2]*n[2] + tau[1]*(u[1]-uhat[1]);
    fb = [0.0 tm tm tm];
    return fb;
end
function initu(x, mu, eta)
    u0 = 0.0;
    return u0;
end
