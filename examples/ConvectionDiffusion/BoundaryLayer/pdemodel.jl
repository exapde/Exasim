function flux(u, q, w, v, x, t, mu, eta)
    f = mu[1]*q + u[1]*mu[2:3];
    return f;
end
function source(u, q, w, v, x, t, mu, eta)
    x1 = x[1];
    x2 = x[2];
    v1 = mu[2];
    v2 = mu[3];
    s = (v1*x2*exp(v1)*exp(v2)*(exp(v2*(x2 - 1)) - 1)*(exp(v1*(x1 - 1)) + v1*x1*exp(v1*(x1 - 1)) - 1))/((exp(v1) - 1)*(exp(v2) - 1)) - (v2*x1*exp(v1 + v2*x2)*(exp(v1*(x1 - 1)) - 1)*(v2*x2 + 2))/((exp(v1) - 1)*(exp(v2) - 1)) - (v1*x2*exp(v2 + v1*x1)*(exp(v2*(x2 - 1)) - 1)*(v1*x1 + 2))/((exp(v1) - 1)*(exp(v2) - 1)) + (v2*x1*exp(v1)*exp(v2)*(exp(v1*(x1 - 1)) - 1)*(exp(v2*(x2 - 1)) + v2*x2*exp(v2*(x2 - 1)) - 1))/((exp(v1) - 1)*(exp(v2) - 1));
    return s;
end
function ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    ub = 0.0;
    return ub;
end
function fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    f = flux(u, q, w, v, x, t, mu, eta);
    fb = f[1]*n[1] + f[2]*n[2] + tau[1]*(u[1]-uhat[1]);
    return fb;
end
function initu(x, mu, eta)
    u0 = 0.0;
    return u0;
end
