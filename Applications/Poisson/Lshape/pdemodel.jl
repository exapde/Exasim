# Poisson equation with homogenous Dirichlet condition on a unit square
# q + \nabla u = 0 in \Omega
# \nabla dot flux = source in \Omega
# u = 0 on \partial Omega
# flux = 2*param*q
# source = 2*pi*pi*sin(pi*x1)*sin(pi*x2);

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
