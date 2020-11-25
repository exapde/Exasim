function pde = pdemodel
pde.mass = @mass;
pde.flux = @flux;
pde.source = @source;
pde.fbou = @fbou;
pde.ubou = @ubou;
pde.initu = @initu;
end

function m = mass(u, q, w, v, x, t, mu, eta)
m = sym([1.0; 1.0; 1.0; 1.0]); 
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
    rE = u(4);
    rx = q(1);
    rux = q(2);
    rvx = q(3);
    rEx = q(4);
    ry = q(5);
    ruy = q(6);
    rvy = q(7);
    rEy = q(8);
    r1 = 1/r;
    uv = ru*r1;
    vv = rv*r1;
    E = rE*r1;
    ke = 0.5*(uv*uv+vv*vv);
    p = gam1*(rE-r*ke);
    h = E+p*r1;
    fi = [ru, ru*uv+p, rv*uv, ru*h, ...
            rv, ru*vv, rv*vv+p, rv*h];
    ux = (rux - rx*uv)*r1;
    vx = (rvx - rx*vv)*r1;
    kex = uv*ux + vv*vx;
    px = gam1*(rEx - rx*ke - r*kex);
    Tx = gam*M2*(px*r - p*rx)*r1^2;
    uy = (ruy - ry*uv)*r1;
    vy = (rvy - ry*vv)*r1;
    key = uv*uy + vv*vy;
    py = gam1*(rEy - ry*ke - r*key);
    Ty = gam*M2*(py*r - p*ry)*r1^2;
    txx = Re1*c23*(2*ux - vy);
    txy = Re1*(uy + vx);
    tyy = Re1*c23*(2*vy - ux);
    fv = [0, txx, txy, uv*txx + vv*txy + fc*Tx, ...
          0, txy, tyy, uv*txy + vv*tyy + fc*Ty];
    f = fi+fv;
    f = reshape(f,[4,2]);    
end

function s = source(u, q, w, v, x, t, mu, eta)
    gravity = mu(5);
    s = [sym(0.0); sym(0.0); -gravity*u(1); -gravity*u(3)];
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)

    f = flux(u, q, w, v, x, t, mu, eta);
    fi = f(:,1)*n(1) + f(:,2)*n(2) + tau*(u-uhat); % numerical flux at freestream boundary
    fw = fi;
    fw(1) = 0.0;   % zero velocity 
    fb = [fw fw];  % wall and freestream boundary conditions    
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
    utw1(2:3) = 0.0;
    utw1(4) = u(1)*Tbot/(gam*gam1*M2);

    utw2 = u;
    utw2(2:3) = 0.0;
    utw2(4) = u(1)*Ttop/(gam*gam1*M2);

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
    z = x(2);
    r = rbot*exp(-gam*M2*gravity*z/Tbot);
    p = pbot*exp(-gam*M2*gravity*z/Tbot);
    u0 = [ r, 0.0, 0.0, p/(gam-1.0)];    
end


