function pde = pdemodel
pde.mass = @mass;
pde.flux = @flux;
pde.source = @source;
pde.fbou = @fbou;
pde.ubou = @ubou;
pde.initu = @initu;
pde.stab = @stab;
end

function m = mass(u, q, w, v, x, t, mu, eta)
m = sym([1.0; 1.0; 1.0; 1.0]); 
end

function f = flux(u, q, w, v, x, t, mu, eta)
    gam = mu(1);
    Re = mu(2);
    Pe = mu(3);
    Minf = mu(4);
    M2 = Minf^2;
    c23 = 2.0/3.0;
    
    r = u(1);
    ruv = u(2);
    rvv = u(3);
    rT = u(4);
    
    rho = exp(r);
    sr = sqrt(rho);
    sr1 = 1/sr;
    
    uv = ruv*sr1;
    vv = rvv*sr1;
    T = rT*sr1;
    
    p = rT/(gam*M2);
    
    rx = q(1);
    rux = q(2);
    rvx = q(3);
    rTx = q(4);
    ry = q(5);
    ruy = q(6);
    rvy = q(7);
    rTy = q(8);
    
    ux = sr1*rux - 0.5*rx*uv;
    uy = sr1*ruy - 0.5*ry*uv;
    vx = sr1*rvx - 0.5*rx*vv;
    vy = sr1*rvy - 0.5*ry*vv;
    Tx = sr1*rTx - 0.5*rx*T;
    Ty = sr1*rTy - 0.5*ry*T;   

    fi = [r*uv, ruv*uv+p, rvv*uv, uv*rT, ...
            r*vv, ruv*vv, rvv*vv+p, vv*rT];
        
    % Viscosity
    mustar = sqrt(T);
    kstar = T^0.75;
    nu = mustar*sr1/Re;
    fc = kstar*sr1*gam/Pe;
    
    txx = nu*c23*(2*ux - vy);
    txy = nu*(uy + vx);
    tyy = nu*c23*(2*vy - ux);
    
    fv = [0, txx, txy, fc*Tx, ...
          0, txy, tyy, fc*Ty];
    f = fi+fv;
    f = reshape(f,[4,2]);    
end

function s = source(u, q, w, v, x, t, mu, eta)
    gam = mu(1);
    gam1 = gam - 1.0;
    Re = mu(2);
    Pe = mu(3);
    Minf = mu(4);
    Re1 = 1/Re;
    M2 = Minf^2;
    c23 = 2.0/3.0;
    
    r = u(1);
    ruv = u(2);
    rvv = u(3);
    rT = u(4);
    
    r_1 = r-1;
    rho = exp(r);
    sr = sqrt(rho);
    sr1 = 1/sr;
    
    uv = ruv*sr1;
    vv = rvv*sr1;
    T = rT*sr1;
    
    p = rT/(gam*M2);
    
    rx = -q(1);
    rux = -q(2);
    rvx = -q(3);
    rTx = -q(4);
    ry = -q(5);
    ruy = -q(6);
    rvy = -q(7);
    rTy = -q(8);
    
    ux = sr1*rux - 0.5*rx*uv;
    uy = sr1*ruy - 0.5*ry*uv;
    vx = sr1*rvx - 0.5*rx*vv;
    vy = sr1*rvy - 0.5*ry*vv;
    Tx = sr1*rTx - 0.5*rx*T;
    Ty = sr1*rTy - 0.5*ry*T; 
    
    % Viscosity
    mustar = sqrt(T);
    kstar = T^0.75;
    nu = mustar*sr1/Re;
    fc = kstar*sr1*gam/Pe;
    
    txx = nu*c23*(2*ux - vy);
    txy = nu*(uy + vx);
    tyy = nu*c23*(2*vy - ux);

    z = sqrt(x(1)^2 + x(2)^2);
    R0 = mu(10);
    gravity0 = mu(5);
    gravity = gravity0*R0^2/(z^2);
    omega = mu(6);
    ax = -gravity*x(1)/z + 2*omega*vv + x(1)*omega^2;
    ay = -gravity*x(2)/z - 2*omega*uv + x(2)*omega^2;
    
    div = ux + vy;
    SigmadV = (txx*ux + txy*uy + txy*vx + tyy*vy)*(gam*gam1*M2);
    Sigmadrx = txx*rx + txy*ry;
    Sigmadry = txy*rx + tyy*ry;
    dTdr = fc*(Tx*rx + Ty*ry);
    s_EUV = EUVsource(u, x, t, mu, eta);
    
    s = [r_1*div; ...
        sr*ax + 0.5*div*ruv - 0.5*p*rx + 0.5*Sigmadrx; ...
        sr*ay + 0.5*div*rvv - 0.5*p*ry + 0.5*Sigmadry; ...
        sr*s_EUV + (3/2-gam)*rT*div + 0.5*dTdr + SigmadV];
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    tau = gettau(uhat, mu, eta, n);

    f = flux(u, q, w, v, x, t, mu, eta);
    fh = f(:,1)*n(1) + f(:,2)*n(2); % numerical flux at freestream boundary
    fw = fh;
    fw(1) = tau(1)*(u(1)-uhat(1));
    fw(2) = fw(2) + tau(2)*(u(2)-uhat(2));
    fw(3) = fw(3) + tau(3)*(u(3)-uhat(3));
    fw(4) = fw(4) + tau(4)*(u(4)-uhat(4));
    
    % Inviscid wall
    gam = mu(1);
    Minf = mu(4);
    M2 = Minf^2;
    
    r = u(1);
    ruv = u(2);
    rvv = u(3);
    rT = u(4);
    
    rho = exp(r);
    sr = sqrt(rho);
    sr1 = 1/sr;
    
    uv = ruv*sr1;
    vv = rvv*sr1;
    
    p = rT/(gam*M2);  

    fi = [r*uv, ruv*uv+p, rvv*uv, uv*rT, ...
            r*vv, ruv*vv, rvv*vv+p, vv*rT];
    fi = reshape(fi,[4,2]);    
    ft = fi(:,1)*n(1) + fi(:,2)*n(2);    
    
    ft(1) = ft(1) + tau(1)*(u(1)-uhat(1));
    ft(2) = ft(2) + tau(2)*(u(2)-uhat(2));
    ft(3) = ft(3) + tau(3)*(u(3)-uhat(3));
    ft(4) = ft(4) + tau(4)*(u(4)-uhat(4));
    
    fb = [fw ft];
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    Tbot = mu(8);

    % Isothermal Wall
    r = u(1);
    rho = exp(r);
    sr = sqrt(rho);
    
    utw1 = u;
    utw1(2:3) = 0.0;
    utw1(4) = sr*Tbot;

    % Inviscid wall
    utw2 = u;
    
    ub = [utw1 utw2];
end


function u0 = initu(x, mu, eta)
    radius = sqrt(x(1)^2 + x(2)^2);
    
    gam = mu(1);
    Minf = mu(4);
    M2 = Minf^2;
    Fr2 = mu(5);
    omega = mu(6);

    Tbot = mu(8);
    Ttop = mu(9);
    R0 = mu(10);
    Ldim = mu(14);
    h0 = 40000/Ldim;

    a0 = gam*M2*(-Fr2 + omega^2*R0);
    
    T = Ttop - (Ttop-Tbot)*exp(-(radius-R0)/h0);
    logp_p0 = a0*h0/Ttop*log(1+Ttop/Tbot*(exp((radius-R0)/h0)-1));
    r = logp_p0 - log(T);
    rho = exp(r);
    srT = sqrt(rho)*T;
    
    u0 = sym([r; 0.0; 0.0; srT]);
end


function ftau = stab(u1, q1, w1, v1, x, t, mu, eta, uhat, n, tau, u2, q2, w2, v2) 
    uhat = 0.5*(u1+u2);
    tau = gettau(uhat, mu, eta, n);
    
    ftau(1) = tau(1)*(u1(1) - u2(1));
    ftau(2) = tau(2)*(u1(2) - u2(2));
    ftau(3) = tau(3)*(u1(3) - u2(3));
    ftau(4) = tau(4)*(u1(4) - u2(4));
end


function   tau = gettau(uhat, mu, eta, n)

    gam = mu(1);
    Re = mu(2);
    Pe = mu(3);
    
    r = uhat(1);
    rT = uhat(4);
    
    rho = exp(r);
    sr = sqrt(rho);
    sr1 = 1/sr;
    T = rT*sr1;
    
    tauA = mu(18);
    
    % Viscosity
    mustar = sqrt(T);
    kstar = T^0.75;
    tauDv = mustar*sr1/Re;
    tauDT = kstar*sr1*gam/Pe;
    
    tau = 0*uhat;

    tau(1) = tauA;
    tau(2) = tauA + tauDv;
    tau(3) = tauA + tauDv;
    tau(4) = tauA + tauDT;

end