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
    ub = [u u];
    return ub;
end
function fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    #ui = sin(mu[2]*x[1]+mu[3]*x[2]-mu[4]*sqrt(mu[2]*mu[2]+mu[3]*mu[3])*t);
    uix = mu[2]*cos(mu[2]*x[1]+mu[3]*x[2]-mu[4]*sqrt(mu[2]*mu[2]+mu[3]*mu[3])*t);
    uiy = mu[3]*cos(mu[2]*x[1]+mu[3]*x[2]-mu[4]*sqrt(mu[2]*mu[2]+mu[3]*mu[3])*t);
    fb = [u uix*n[1]+uiy*n[2]]; # Neumann conditions
    return fb;
end
function initu(x, mu, eta)
    u0 = 0.0;
    return u0;
end
function initq(x, mu, eta)
    q0 = [0.0 0.0];
end
function initw(x, mu, eta)
    w0 = 0.0;
    return w0;
end
