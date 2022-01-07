function mass(u, q, w, v, x, t, mu, eta)
    m = [1.0, 1.0, 1.0, 1.0, 1.0];
    return m;
end
function flux(u, q, w, v, x, t, mu, eta)
    gam = mu[1];    
    gam1 = gam - 1.0;
    Re = mu[2];
    Pr = mu[3];
    Minf = mu[4];
    Re1 = 1/Re;
    M2 = Minf^2;
    c23 = 2.0/3.0;
    fc = 1/(gam1*M2*Re*Pr);
    r = u[1];
    ru = u[2];
    rv = u[3];
    rw = u[4];
    rE = u[5];
    rx = q[1];
    rux = q[2];
    rvx = q[3];
    rwx = q[4];
    rEx = q[5];
    ry = q[6];
    ruy = q[7];
    rvy = q[8];
    rwy = q[9];
    rEy = q[10];
    rz = q[11];
    ruz = q[12];
    rvz = q[13];
    rwz = q[14];
    rEz = q[15];
    r1 = 1/r;
    uv = ru*r1;
    vv = rv*r1;
    wv = rw*r1;
    E = rE*r1;
    ke = 0.5*(uv*uv+vv*vv+wv*wv);
    p = gam1*(rE-r*ke);
    h = E+p*r1;
    fi = [ru, ru*uv+p, rv*uv, rw*uv, ru*h, rv, ru*vv, rv*vv+p, rw*vv, rv*h, rw, ru*wv, rv*wv, rw*wv+p, rw*h];
    ux = (rux - rx*uv)*r1;
    vx = (rvx - rx*vv)*r1;
    wx = (rwx - rx*wv)*r1;
    kex = uv*ux + vv*vx + wv*wx;
    px = gam1*(rEx - rx*ke - r*kex);
    Tx = gam*M2*(px*r - p*rx)*r1^2;
    uy = (ruy - ry*uv)*r1;
    vy = (rvy - ry*vv)*r1;
    wy = (rwy - ry*wv)*r1;
    key = uv*uy + vv*vy + wv*wy;
    py = gam1*(rEy - ry*ke - r*key);
    Ty = gam*M2*(py*r - p*ry)*r1^2;
    uz = (ruz - rz*uv)*r1;
    vz = (rvz - rz*vv)*r1;
    wz = (rwz - rz*wv)*r1;
    kez = uv*uz + vv*vz + wv*wz;
    pz = gam1*(rEz - rz*ke - r*kez);
    Tz = gam*M2*(pz*r - p*rz)*r1^2;
    txx = Re1*c23*(2*ux - vy - wz);
    txy = Re1*(uy + vx);
    txz = Re1*(uz + wx);
    tyy = Re1*c23*(2*vy - ux - wz);
    tyz = Re1*(vz + wy);
    tzz = Re1*c23*(2*wz - ux - vy);
    fv = [0, txx, txy, txz, uv*txx + vv*txy + wv*txz + fc*Tx, 0, txy, tyy, tyz, uv*txy + vv*tyy + wv*tyz + fc*Ty, 0, txz, tyz, tzz, uv*txz + vv*tyz + wv*tzz + fc*Tz];
    f = fi+fv;
    f = reshape(f,(5,3));    
    return f;
end
function source(u, q, w, v, x, t, mu, eta)
    s = [0.0, 0.0, 0.0, 0.0, 0.0];
    return s;
end
function ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    ub = [0.0, 0.0, 0.0, 0.0, 0.0];
    return ub;
end
function fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    f = flux(u, q, w, v, x, t, mu, eta);
    fb = f[:,1]*n[1] + f[:,2]*n[2] + f[:,3]*n[3] + tau[1]*(u-uhat);
    return fb;
end
function initu(x, mu, eta)
    L = 2*pi;
    x1 = x[1];
    x2 = x[2];
    x3 = x[3];
    gam = mu[1];
    Minf = mu[4];
    p0 = 1/(gam*Minf^2);
    
    u01 = 1.0;
    u02 = sin(2*pi*(x1-L/2)/L) * cos(2*pi*(x2-L/2)/L) * cos(2*pi*(x3-L/2)/L);
    u03 = - cos(2*pi*(x1-L/2)/L) * sin(2*pi*(x2-L/2)/L) * cos(2*pi*(x3-L/2)/L);
    u04 = 0.0;
    p = p0 + (1/16) * (cos(2*2*pi*(x1-L/2)/L) + cos(2*2*pi*(x2-L/2)/L)) * (cos(2*2*pi*(x3-L/2)/L) + 2.0);
    u05 = 0.5*(u02^2+u03^2+u04^2) / u01 + p / (gam-1);

    u0 = [u01, u02, u03, u04, u05];
    return u0;
end

