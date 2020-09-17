# Poisson equation with homogenous Dirichlet condition on a unit square
# q + \nabla u = 0 in \Omega
# \nabla dot flux = source in \Omega
# u = 0 on \partial Omega
# flux = 2*param*q
# source = 2*pi*pi*sin(pi*x)*sin(pi*y);
# udg = (u, qx, qy) with udg[1] = u, udg[2] = qx, udg[3] = qy

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
