function mass(u, q, w, v, x, t, mu, eta)
    m = [1.0, 1.0, 1.0, 1.0];
    return m;
end
function flux(u, q, w, v, x, t, mu, eta)
    gam = mu[1];
    gam1 = gam - 1.0;
    r = u[1];
    ru = u[2];
    rv = u[3];
    rE = u[4];
    r1 = 1/r;
    uv = ru*r1;
    vv = rv*r1;
    E = rE*r1;
    p = gam1*(rE-r*0.5*(uv*uv+vv*vv));
    h = E+p*r1;
    f = [ru, ru*uv+p, rv*uv, ru*h, rv, ru*vv, rv*vv+p, rv*h];
    f = reshape(f,(4,2));
    return f;
end
function source(u, q, w, v, x, t, mu, eta)
    s = [0.0*u[1], 0.0*u[1], 0.0*u[1], 0.0*u[1]];
    return s;
end
function ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    ub = [0.0, 0.0, 0.0, 0.0];
    return ub;
end
function fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    f = flux(u, q, w, v, x, t, mu, eta);
    fb = f[:,1]*n[1] + f[:,2]*n[2] + tau[1]*(u-uhat);
    return fb;
end
function initu(x, mu, eta)

    t = 0.0;
    phi = 5.0;
    x1 = x[1];
    x2 = x[2];
    gam = mu[1];
    M_ref = mu[2];

    r = sqrt((x1-t)^2 + x2^2);
    u01 = (1 - ((gam-1)/(16*pi^2)) * phi^2 * exp(2*(1-r^2)) )^(1/(gam-1));
    u02 = u01 * (1 - M_ref^(-1) * phi * exp(1-r^2) * x2/(2*pi));
    u03 = u01 * M_ref^(-1) * phi * exp(1-r^2) * x1/(2*pi);
    p = u01^gam/(gam*M_ref^2);
    u04 = p/(gam-1) + 0.5 * (u02*u02 + u03*u03) / u01;

    u0 = [u01, u02, u03, u04];
    return u0;
end

function fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau)
fb0 = tau[1]*(0*u-uhat);
fb = [fb0 fb0 fb0 fb0];
return fb
end
