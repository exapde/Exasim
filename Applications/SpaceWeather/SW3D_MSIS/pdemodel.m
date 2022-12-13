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
    m = sym([1.0; 1.0; 1.0; 1.0; 1.0]); 
end

function f = flux(u, q, w, v, x, t, mu, eta)
    fi = fluxInviscid(u,x,mu,eta);
    fv = fluxViscous(u,q,x,mu,eta);
    
    f = fi+fv;
end

function f = fluxWall(u, q, w, v, x, t, mu, eta)
    fi = fluxInviscidWall(u,x,mu,eta);
    fv = fluxViscousWall(u,q,x,mu,eta);
    
    f = fi+fv;
end

function s = source(u, q, w, v, x, t, mu, eta)
    gam = mu(1);
    gam1 = gam - 1.0;
    Gr = mu(2);
    Pr = mu(3);
    Fr = mu(4);
    c23 = 2.0/3.0;
    
    [mw,dmdr] = weightedMass3D(x,mu,eta);
    
    r = u(1);
    ruv = u(2);
    rvv = u(3);
    rwv = u(4);
    rT = u(5);
    
    r_1 = r-1;
    rho = exp(r);
    sr = sqrt(rho);
    sr1 = 1/sr;
    
    uv = ruv*sr1;
    vv = rvv*sr1;
    wv = rwv*sr1;
    T = rT*sr1;
    
    Tmin = 0.1;
    alpha = 1e3;
    T = Tmin + lmax(T-Tmin,alpha);
    rT = sr*T;
    
    p = rT/(gam*mw);
    
    rx = -q(1);
    rux = -q(2);
    rvx = -q(3);
    rwx = -q(4);
    rTx = -q(5);
    ry = -q(6);
    ruy = -q(7);
    rvy = -q(8);
    rwy = -q(9);
    rTy = -q(10);
    rz = -q(11);
    ruz = -q(12);
    rvz = -q(13);
    rwz = -q(14);
    rTz = -q(15);
    
    ux = sr1*rux - 0.5*rx*uv;
    uy = sr1*ruy - 0.5*ry*uv;
    uz = sr1*ruz - 0.5*rz*uv;
    vx = sr1*rvx - 0.5*rx*vv;
    vy = sr1*rvy - 0.5*ry*vv;
    vz = sr1*rvz - 0.5*rz*vv;
    wx = sr1*rwx - 0.5*rx*wv;
    wy = sr1*rwy - 0.5*ry*wv;
    wz = sr1*rwz - 0.5*rz*wv;
    Tx = sr1*rTx - 0.5*rx*T;
    Ty = sr1*rTy - 0.5*ry*T;
    Tz = sr1*rTz - 0.5*rz*T;
    
    % Temperature sensor
    dT = atan(alpha*(T - Tmin))/pi + (alpha*(T - Tmin))/(pi*(alpha^2*(T - Tmin)^2 + 1)) + 1/2;
    Tx = Tx*dT;
    Ty = Ty*dT;
    Tz = Tz*dT;
    
    % Viscosity
    expMu = mu(12);
    expKappa = mu(13);
    nuEddy = mu(14);
    alphaEddy = mu(15);
    
    mustar = T^expMu;
    k0 = ThermalConductivity3D(x,mu,eta);
    kstar = k0*T^expKappa;
    nu = (mustar*sr1 + sr*nuEddy)/sqrt(gam*Gr);
    fc = (kstar*sr1 + sr*alphaEddy)*mw*sqrt(gam/Gr)/Pr;
    
    txx = nu*c23*(2*ux - vy - wz);
    tyy = nu*c23*(2*vy - ux - wz);
    tzz = nu*c23*(2*wz - ux - vy);
    txy = nu*(uy + vx);
    txz = nu*(uz + wx);
    tyz = nu*(vz + wy);

    z = sqrt(x(1)^2 + x(2)^2 + x(3)^2);
    R0 = mu(16);
    gravity0 = 1/gam;
    gravity = gravity0*R0^2/(z^2);
    omega = Fr/sqrt(gam);
    
    ax = -gravity*x(1)/z + 2*omega*vv + x(1)*omega^2;
    ay = -gravity*x(2)/z - 2*omega*uv + x(2)*omega^2;
    az = -gravity*x(3)/z;
    
    div = ux + vy + wz;
    SigmadV = (txx*ux + txy*(uy+vx) + txz*(uz+wx) + tyy*vy + tyz*(vz+wy) + tzz*wz)*mw*gam*gam1;
    Sigmadrx = txx*rx + txy*ry + txz*rz;
    Sigmadry = txy*rx + tyy*ry + tyz*rz;
    Sigmadrz = txz*rx + tyz*ry + tzz*rz;
    dTdr = fc*(Tx*rx + Ty*ry + Tz*rz);
    s_EUV = EUVsource3D(u, x, t, mu, eta);
    
    %source due to mass variation;
    qcv = (dmdr/mw)*((rT*vx - fc*Tx)*x(1) + (rT*vy - fc*Ty)*x(2) + (rT*vz - fc*Tz)*x(3))/z;
      
    s = [r_1*div; ...
        sr*ax + 0.5*div*ruv - 0.5*p*rx + 0.5*Sigmadrx; ...
        sr*ay + 0.5*div*rvv - 0.5*p*ry + 0.5*Sigmadry; ...
        sr*az + 0.5*div*rwv - 0.5*p*rz + 0.5*Sigmadrz; ...
        sr*s_EUV + (3/2-gam)*rT*div + 0.5*dTdr + SigmadV + qcv];
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    tau = gettau(uhat, mu, eta, x, n);

    f = fluxWall(u, q, w, v, x, t, mu, eta);
    fh = f(:,1)*n(1) + f(:,2)*n(2) + f(:,3)*n(3); % numerical flux at freestream boundary
    fw = fh;
    fw(1) = fw(1) + tau(1)*(u(1)-uhat(1));
    fw(2) = fw(2) + tau(2)*(u(2)-uhat(2));
    fw(3) = fw(3) + tau(3)*(u(3)-uhat(3));
    fw(4) = fw(4) + tau(4)*(u(4)-uhat(4));
    fw(5) = fw(5) + tau(5)*(u(5)-uhat(5));
    
    % Inviscid wall
    fi = fluxInviscid(u,x,mu,eta);  
    ft = fi(:,1)*n(1) + fi(:,2)*n(2) + fi(:,3)*n(3);    
    
    ft(1) = ft(1) + tau(1)*(u(1)-uhat(1));
    ft(2) = ft(2) + tau(2)*(u(2)-uhat(2));
    ft(3) = ft(3) + tau(3)*(u(3)-uhat(3));
    ft(4) = ft(4) + tau(4)*(u(4)-uhat(4));
    ft(5) = ft(5) + tau(5)*(u(5)-uhat(5));
    
    fb = [fw ft];
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    Tbot = 1.0;

    % Isothermal Wall
    r = u(1);
    rho = exp(r);
    sr = sqrt(rho);
    
    utw1 = u;
    utw1(2:4) = 0.0;
    utw1(5) = sr*Tbot;

    % Inviscid wall
    utw2 = u;
    
    ub = [utw1 utw2];
end

function u0 = initu(x, mu, eta)
    radius = sqrt(x(1)^2 + x(2)^2 + x(3)^2);    
    Fr = mu(4);
    
    mw = weightedMass3D(x,mu,eta);

    Tbot = 1.0;
    Ttop = 5.0;
    R0 = mu(16);
    H = mu(18);
    h0 = 35000/H;

    a0 = -1 + Fr^2*R0;
    
    T = Ttop - (Ttop-Tbot)*exp(-(radius-R0)/h0);
    logp_p0 = a0*mw*h0/Ttop*log(1+Ttop/Tbot*(exp((radius-R0)/h0)-1));
    r = logp_p0 - log(T) + log(mw);
    rho = exp(r);
    srT = sqrt(rho)*T;
    
    u0 = sym([r; 0.0; 0.0; 0.0; srT]);
end

function ftau = stab(u1, q1, w1, v1, x, t, mu, eta, uhat, n, tau, u2, q2, w2, v2) 
    uhat = 0.5*(u1+u2);
    tau = gettau(uhat, mu, eta, x, n);
    
    ftau(1) = tau(1)*(u1(1) - u2(1));
    ftau(2) = tau(2)*(u1(2) - u2(2));
    ftau(3) = tau(3)*(u1(3) - u2(3));
    ftau(4) = tau(4)*(u1(4) - u2(4));
    ftau(5) = tau(5)*(u1(5) - u2(5));
end

function fi = fluxInviscid(u,x,mu,eta)
    gam = mu(1);
    
    r = u(1);
    ruv = u(2);
    rvv = u(3);
    rwv = u(4);
    rT = u(5);
    
    rho = exp(r);
    sr = sqrt(rho);
    sr1 = 1/sr;
    
    uv = ruv*sr1;
    vv = rvv*sr1;
    wv = rwv*sr1;
    T = rT*sr1;    
    
    Tmin = 0.1;
    alpha = 1e3;
    T = Tmin + lmax(T-Tmin,alpha);
    rT = sr*T;
    
    mw = weightedMass3D(x,mu,eta);
    p = rT/(gam*mw);
    
    fi = [r*uv, ruv*uv+p, rvv*uv, rwv*uv, uv*rT, ...
      r*vv, ruv*vv, rvv*vv+p, rwv*vv, vv*rT, ...
      r*wv, ruv*wv, rvv*wv, rwv*wv+p, wv*rT];

    fi = reshape(fi,[5,3]); 
end

function fi = fluxInviscidWall(u,x,mu,eta)
    gam = mu(1);
    
    r = u(1);
    rho = exp(r);
    sr = sqrt(rho);
    
    Tbot = 1.0;    
    rT = sr*Tbot;
  
    mw = weightedMass3D(x,mu,eta);
    p = rT/(gam*mw);
    
    fi = [0.0, p, 0.0, 0.0, 0.0, ...
      0.0, 0.0, p, 0.0, 0.0, ...
      0.0, 0.0, 0.0, p, 0.0];

    fi = reshape(fi,[5,3]); 
end

function fv = fluxViscous(u, q, x, mu, eta)
    gam = mu(1);
    Gr = mu(2);
    Pr = mu(3);
    c23 = 2.0/3.0;
    mw = weightedMass3D(x,mu,eta);
    
    r = u(1);
    ruv = u(2);
    rvv = u(3);
    rwv = u(4);
    rT = u(5);
    
    rho = exp(r);
    sr = sqrt(rho);
    sr1 = 1/sr;
    
    uv = ruv*sr1;
    vv = rvv*sr1;
    wv = rwv*sr1;
    T = rT*sr1;    
    
    Tmin = 0.1;
    alpha = 1e3;
    T = Tmin + lmax(T-Tmin,alpha);
      
    rx  = q(1);
    rux = q(2);
    rvx = q(3);
    rwx = q(4);
    rTx = q(5);
    ry  = q(6);
    ruy = q(7);
    rvy = q(8);
    rwy = q(9);
    rTy = q(10);
    rz  = q(11);
    ruz = q(12);
    rvz = q(13);
    rwz = q(14);
    rTz = q(15);
    
    ux = sr1*rux - 0.5*rx*uv;
    uy = sr1*ruy - 0.5*ry*uv;
    uz = sr1*ruz - 0.5*rz*uv;
    vx = sr1*rvx - 0.5*rx*vv;
    vy = sr1*rvy - 0.5*ry*vv;
    vz = sr1*rvz - 0.5*rz*vv;
    wx = sr1*rwx - 0.5*rx*wv;
    wy = sr1*rwy - 0.5*ry*wv;
    wz = sr1*rwz - 0.5*rz*wv;
    Tx = sr1*rTx - 0.5*rx*T;
    Ty = sr1*rTy - 0.5*ry*T;
    Tz = sr1*rTz - 0.5*rz*T;
    
    % Temperature sensor
    dT = atan(alpha*(T - Tmin))/pi + (alpha*(T - Tmin))/(pi*(alpha^2*(T - Tmin)^2 + 1)) + 1/2;
    Tx = Tx*dT;
    Ty = Ty*dT;
    Tz = Tz*dT;
        
    % Viscosity
    expMu = mu(12);
    expKappa = mu(13);
    nuEddy = mu(14);
    alphaEddy = mu(15);
    
    mustar = T^expMu;
    k0 = ThermalConductivity3D(x,mu,eta);
    kstar = k0*T^expKappa;
    nu = (mustar*sr1 + sr*nuEddy)/sqrt(gam*Gr);
    fc = (kstar*sr1 + sr*alphaEddy)*mw*sqrt(gam/Gr)/Pr;
    
    txx = nu*c23*(2*ux - vy - wz);
    tyy = nu*c23*(2*vy - ux - wz);
    tzz = nu*c23*(2*wz - ux - vy);
    txy = nu*(uy + vx);
    txz = nu*(uz + wx);
    tyz = nu*(vz + wy);
    
    fv = [0, txx, txy, txz, fc*Tx, ...
          0, txy, tyy, tyz, fc*Ty, ...
          0, txz, tyz, tzz, fc*Tz];
    fv = reshape(fv,[5,3]);    
end

function fv = fluxViscousWall(u, q, x, mu, eta)
    gam = mu(1);
    Gr = mu(2);
    Pr = mu(3);
    c23 = 2.0/3.0;
    mw = weightedMass3D(x,mu,eta);
    
    r = u(1);
    
    rho = exp(r);
    sr = sqrt(rho);
    sr1 = 1/sr;
    
    uv = 0.0;
    vv = 0.0;
    wv = 0.0;
    Tbot = 1.0;    
      
    rx  = q(1);
    rux = q(2);
    rvx = q(3);
    rwx = q(4);
    rTx = q(5);
    ry  = q(6);
    ruy = q(7);
    rvy = q(8);
    rwy = q(9);
    rTy = q(10);
    rz  = q(11);
    ruz = q(12);
    rvz = q(13);
    rwz = q(14);
    rTz = q(15);
    
    ux = sr1*rux - 0.5*rx*uv;
    uy = sr1*ruy - 0.5*ry*uv;
    uz = sr1*ruz - 0.5*rz*uv;
    vx = sr1*rvx - 0.5*rx*vv;
    vy = sr1*rvy - 0.5*ry*vv;
    vz = sr1*rvz - 0.5*rz*vv;
    wx = sr1*rwx - 0.5*rx*wv;
    wy = sr1*rwy - 0.5*ry*wv;
    wz = sr1*rwz - 0.5*rz*wv;
    Tx = sr1*rTx - 0.5*rx*Tbot;
    Ty = sr1*rTy - 0.5*ry*Tbot;
    Tz = sr1*rTz - 0.5*rz*Tbot;
        
    % Viscosity
    expMu = mu(12);
    expKappa = mu(13);
    nuEddy = mu(14);
    alphaEddy = mu(15);
    
    mustar = Tbot^expMu;
    k0 = ThermalConductivity3D(x,mu,eta);
    kstar = k0*Tbot^expKappa;
    nu = (mustar*sr1 + sr*nuEddy)/sqrt(gam*Gr);
    fc = (kstar*sr1 + sr*alphaEddy)*mw*sqrt(gam/Gr)/Pr;
    
    txx = nu*c23*(2*ux - vy - wz);
    tyy = nu*c23*(2*vy - ux - wz);
    tzz = nu*c23*(2*wz - ux - vy);
    txy = nu*(uy + vx);
    txz = nu*(uz + wx);
    tyz = nu*(vz + wy);
    
    fv = [0, txx, txy, txz, fc*Tx, ...
          0, txy, tyy, tyz, fc*Ty, ...
          0, txz, tyz, tzz, fc*Tz];
    fv = reshape(fv,[5,3]);    
end

function   tau = gettau(uhat, mu, eta, x, n)

    gam = mu(1);
    Gr = mu(2);
    Pr = mu(3);
    mw = weightedMass3D(x,mu,eta);
    
    r = uhat(1);
    ruv = uhat(2);
    rvv = uhat(3);
    rwv = uhat(4);
    rT = uhat(5);
    
    rho = exp(r);
    sr = sqrt(rho);
    sr1 = 1/sr;
    
    uv = ruv*sr1;
    vv = rvv*sr1;
    wv = rwv*sr1;
    T = rT*sr1;
    
    Tmin = 0.1;
    alpha = 1e3;
    T = Tmin + lmax(T-Tmin,alpha);

%     vn = uv*n(1) + vv*n(2) + wv*n(3);
%     c = sqrt(T);
%     tauA = sqrt(vn*vn) + c;
    tauA = mu(22);

    
    % Viscosity
    expMu = mu(12);
    expKappa = mu(13);
    nuEddy = mu(14);
    alphaEddy = mu(15);
    
    mustar = T^expMu;
    k0 = ThermalConductivity3D(x,mu,eta);
    kstar = k0*T^expKappa;
    tauDv = (mustar*sr1 + sr*nuEddy)/sqrt(gam*Gr);
    tauDT = (kstar*sr1 + sr*alphaEddy)*mw*sqrt(gam/Gr)/Pr;

    tau(1) = tauA;
    tau(2) = tauA + tauDv;
    tau(3) = tauA + tauDv;
    tau(4) = tauA + tauDv;
    tau(5) = tauA + tauDT;
end