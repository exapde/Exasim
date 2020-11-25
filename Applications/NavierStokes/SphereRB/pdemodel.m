function pde = pdemodel
pde.mass = @mass;
pde.flux = @flux;
pde.source = @source;
pde.fbou = @fbou;
pde.ubou = @ubou;
pde.initu = @initu;
end

function m = mass(u, q, w, v, x, t, mu, eta)
m = sym([1.0; 1.0; 1.0; 1.0; 1.0]); 
end

function f = flux(u, q, w, v, x, t, mu, eta)
    gam = mu(1);    
    gam1 = gam - 1.0;
    Re = mu(2);
    Pr = mu(3);
    Minf = mu(4);
    Re1 = 1/Re;
    M2 = Minf^2;
    c23 = 2.0/3.0;
    fc = 1/(gam1*M2*Re*Pr);
    r = u(1);
    ru = u(2);
    rv = u(3);
    rw = u(4);
    rE = u(5);
    rx = q(1);
    rux = q(2);
    rvx = q(3);
    rwx = q(4);
    rEx = q(5);
    ry = q(6);
    ruy = q(7);
    rvy = q(8);
    rwy = q(9);
    rEy = q(10);
    rz = q(11);
    ruz = q(12);
    rvz = q(13);
    rwz = q(14);
    rEz = q(15);
    r1 = 1/r;
    uv = ru*r1;
    vv = rv*r1;
    wv = rw*r1;
    E = rE*r1;
    ke = 0.5*(uv*uv+vv*vv+wv*wv);
    p = gam1*(rE-r*ke);
    h = E+p*r1;
    fi = [ru, ru*uv+p, rv*uv, rw*uv, ru*h, ...
          rv, ru*vv, rv*vv+p, rw*vv, rv*h, ...
          rw, ru*wv, rv*wv, rw*wv+p, rw*h];
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
    fv = [0, txx, txy, txz, uv*txx + vv*txy + wv*txz + fc*Tx, ...
          0, txy, tyy, tyz, uv*txy + vv*tyy + wv*tyz + fc*Ty, ...
          0, txz, tyz, tzz, uv*txz + vv*tyz + wv*tzz + fc*Tz];
    f = fi+fv;
    f = reshape(f,[5,3]);    
end

function s = source(u, q, w, v, x, t, mu, eta)
    r = sqrt(x(1)^2 + x(2)^2 + x(3)^2);
    sca = (x(1)*u(2) + x(2)*u(3) + x(3)*u(4))/r;
    gravity = mu(5);
    s = [sym(0.0); -gravity*u(1)*x(1)/r; ...
                   -gravity*u(1)*x(2)/r; ...
                   -gravity*u(1)*x(3)/r; ...
                   -gravity*sca];
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)

    f = flux(u, q, w, v, x, t, mu, eta);
    fi = f(:,1)*n(1) + f(:,2)*n(2) + f(:,3)*n(3) + tau*(u-uhat); % numerical flux at freestream boundary
    fw = fi;
    fw(1) = 0.0;   % zero velocity 
    fb = [fw fw];  % isothermal wall boundary conditions    
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)

    gam = mu(1);
    gam1 = gam-1.0;
    Minf = mu(4);
    M2 = Minf^2;
    gravity = mu(5);
    rbot = mu(6);
    pbot = mu(7);
    Tbot = mu(8);
    Ttop = mu(9);

    % Isothermal Wall?
    utw1 = u;
    utw1(2:4) = 0.0;
    utw1(5) = u(1)*Tbot/(gam*gam1*M2);

    utw2 = u;
    utw2(2:4) = 0.0;
    utw2(5) = u(1)*Ttop/(gam*gam1*M2);

    ub = [utw1 utw2]; % wall and freestream boundary conditions    
end


function u0 = initu(x, mu, eta)
    gam = mu(1);
    Minf = mu(4);
    M2 = Minf^2;
    gravity = mu(5);
    rbot = mu(6);
    pbot = mu(7);
    Tbot = mu(8);
    R0 = mu(10);
    z = sqrt(x(1)^2 + x(2)^2 + x(3)^2)-R0;
    r = rbot*exp(-gam*M2*gravity*z/Tbot);
    p = pbot*exp(-gam*M2*gravity*z/Tbot);
    u0 = [ r, 0.0, 0.0, 0.0, p/(gam-1.0)];    
end
