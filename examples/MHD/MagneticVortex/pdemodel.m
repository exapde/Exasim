
function pde = pdemodel
pde.mass = @mass;
pde.flux = @flux;
pde.source = @source;
pde.fbou = @fbou;
pde.fbouhdg = @fbouhdg;
pde.ubou = @ubou;
pde.initu = @initu;
end

function m = mass(u, q, w, v, x, t, mu, eta)
m = sym([1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0]); 
end

function f = flux(u, q, w, v, x, t, mu, eta)
    gam = mu(1);
    gam1 = gam - 1.0;    
    r = u(1);
    ru = u(2);
    rv = u(3);
    rE = u(4);
    bx = u(5);
    by = u(6);
    phi = u(7);
    r1 = 1/r;
    uv = ru*r1;
    vv = rv*r1;
    E = rE*r1;
    q = 0.5*(uv*uv+vv*vv);
    b = (bx*bx+by*by)*0.5;
    p = gam1*(rE-r*q-b-0.5*phi^2);
    h = E+p*r1+b*r1;
    uvb = uv*bx+vv*by;
    f = [ru, ru*uv+p+b-bx*bx, rv*uv-by*bx, ru*h-uvb*bx+phi*bx, phi, -(vv*bx-by*uv), bx, ...
            rv, ru*vv-bx*by, rv*vv+p+b-by*by, rv*h-uvb*by+phi*by, -(uv*by-bx*vv), phi, by];
    f = reshape(f,[7,2]);    
end

function s = source(u, q, w, v, x, t, mu, eta)
    s = [sym(0.0); sym(0.0); sym(0.0); sym(0.0); sym(0.0); sym(0.0); sym(0.0)];
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    f = flux(u, q, w, v, x, t, mu, eta);
    fb = f(:,1)*n(1) + f(:,2)*n(2) + tau*(u-uhat);
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    ub = sym([0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0]); 
end

function fb = fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    f = flux(u, q, w, v, x, t, mu, eta);
    fb = f(:,1)*n(1) + f(:,2)*n(2) + tau*(u-uhat);
end

function u0 = initu(x, mu, eta)

    gam = mu(1);
    x1 = x(1);
    x2 = x(2);
    r = sqrt(x1*x1+x2*x2);
    k = 1/(2*pi);    
    u01 = 1;
    u02 = u01-k*x2*exp(0.5*(1-r^2));
    u03 = u01+k*x1*exp(0.5*(1-r^2));
    u05 = -k*x2*exp(0.5*(1-r^2));
    u06 = k*x1*exp(0.5*(1-r^2));
    p = 1 + 1/(4*0.5)*(k^2*(1-2*0.5*r^2)-k^2*u01)*exp(2*0.5*(1-r^2));
    b = u06*u06 + u05*u05;
    v = u02*u02 + u03*u03;
    u04 = p/(gam-1) + u01/2*v + 0.5*b;
    u07 = 0;
    
    u0 = [u01; u02; u03; u04; u05; u06; u07];
end


