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
    m = sym([1.0; 1.0; 1.0]); 
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
    x1 = x(1);

    gam = mu(1);
    gam1 = gam-1;
    Gr = mu(2);
    Pr = mu(3);
    c23 = 2.0/3.0;
    
    [mw,dmdr] = weightedMass(x,mu,eta);
    
    r = u(1);
    srvx = u(2);
    srT = u(3);
    
    rho = exp(r);
    sr = sqrt(rho);
    sr1 = 1/sr;
    r_1 = r-1;
    
    vx = srvx*sr1;
    T = srT*sr1;
    p = srT/(gam*mw);
    
    drdx  = -q(1);
    drvxdx = -q(2);
    drTdx = -q(3);
    
    dvxdx = sr1*drvxdx - 0.5*drdx*vx;
    dTdx = sr1*drTdx - 0.5*drdx*T;
    
    % Viscosity
    expmu = mu(12);
    expkappa = mu(13);
    nuEddy = mu(14);
    alphaEddy = mu(15);
    
    mustar = T^expmu;
%     k0 = ThermalConductivity(x,mu,eta);
    k0 = 1;
    kstar = k0*T^expkappa;
    nu = (mustar*sr1 + sr*nuEddy)/sqrt(gam*Gr);
    fc = (kstar*sr1 + sr*alphaEddy)*mw*sqrt(gam/Gr)/Pr;
    
    trr = nu*c23*2*dvxdx - 2*c23*vx/x1;
    trd = 4*nu*(dvxdx-vx/x1)/x1;
    
    R0 = mu(16);
    gravity0 = 1/gam;
    gravity = gravity0*R0^2/(x1^2);
    Fr = mu(4);
    ar = -gravity + x1*Fr^2/gam;
    
    trp = 2*c23*nu*mw*(dvxdx^2 - 2*vx*dvxdx/x1 + vx^2/x1^2);
    SigmadV = gam*gam1*trp;
    
    q_EUV = EUVsource1D(u, x, t, mu, eta);
    
    %source due to mass variation;
    qcv = (srT*vx - fc*dTdx)*dmdr/mw;
    
    s = [r_1*dvxdx - 2*vx/x1; ...
        sr*ar + 0.5*(dvxdx-2*vx/x1)*srvx - 0.5*p*drdx + 0.5*trr*drdx + 0.5*trd; ...
        sr*q_EUV + (3/2-gam)*srT*dvxdx + 2*(1/2-gam)*srT*vx/x1 + fc*dTdx*(2/x1 + 0.5*drdx) + SigmadV + qcv];
end

function fb = fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    tau = gettau(uhat, mu, eta, x, n);

    f = fluxWall(u, q, w, v, x, t, mu, eta);
    fw = f*n; % numerical flux at freestream boundary
    fw(1) = fw(1) + tau(1)*(u(1)-uhat(1));
    fw(2) = fw(2) + tau(2)*(u(2)-uhat(2));
    fw(3) = fw(3) + tau(3)*(u(3)-uhat(3));
    
    % Inviscid outer boundary
    fi2 = fluxInviscid(u,x,mu,eta);
    ft = fi2*n;    
    ft(1) = ft(1) + tau(1)*(u(1)-uhat(1));
    ft(2) = ft(2) + tau(2)*(u(2)-uhat(2));
    ft(3) = ft(3) + tau(3)*(u(3)-uhat(3));
    
    fb = [fw ft];
end

function ub = ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau)
    Tbot = 1.0;

    % Isothermal Wall
    r = u(1);
    rho = exp(r);
    sr = sqrt(rho);
    
    utw1 = u;
    utw1(2) = 0.0;
    utw1(3) = sr*Tbot;

    % Inviscid wall
    utw2 = u;
        
    ub = [utw1 utw2];
end


function u0 = initu(x, mu, eta)
    x1 = x(1);
    
    Fr = mu(4);

    mw = weightedMass(x,mu,eta);
    
    Tbot = 1.0;
    Ttop = 6.0;
    R0 = mu(16);
    Ldim = mu(18);
    h0 = 35000/Ldim;

    a0 = (-1 + Fr^2*R0);
    
    T = Ttop - (Ttop-Tbot)*exp(-(x1-R0)/h0);
    logp_p0 = a0*mw*h0/Ttop*log(1+Ttop/Tbot*(exp((x1-R0)/h0)-1));
    rtilde = logp_p0 - log(T) + log(mw);
    rho = exp(rtilde);
    srT = sqrt(rho)*T;
    
    u0 = sym([rtilde; 0.0; srT]);
end


function ftau = stab(u1, q1, w1, v1, x, t, mu, eta, uhat, n, tau, u2, q2, w2, v2) 
    uhat = 0.5*(u1+u2);
    tau = gettau(uhat, mu, eta, x, n);
    
    ftau(1) = tau(1)*(u1(1) - u2(1));
    ftau(2) = tau(2)*(u1(2) - u2(2));
    ftau(3) = tau(3)*(u1(3) - u2(3));
end


function fi = fluxInviscid(u,x,mu,eta)
    gam = mu(1);    
    r = u(1);
    srvx = u(2);
    srT = u(3);
    
    rho = exp(r);
    sr = sqrt(rho);
    sr1 = 1/sr;
    
    vx = srvx*sr1;
    
    mw = weightedMass(x,mu,eta);
    p = srT/(gam*mw);

    fi = [r*vx; srvx*vx+p; srT*vx];
end

function fi = fluxInviscidWall(u,x,mu,eta)
    gam = mu(1);
    
    r = u(1);
    rho = exp(r);
    sr = sqrt(rho);
    
    Tbot = 1.0;
    mw = weightedMass(x,mu,eta);
    p = sr*Tbot/(gam*mw);

    fi = [0; p; 0];
end


function fv = fluxViscous(u,q,x,mu,eta)
    x1 = x(1);
    
    gam = mu(1);
    Gr = mu(2);
    Pr = mu(3);
    c23 = 2.0/3.0;
    mw = weightedMass(x,mu,eta);
    
    r = u(1);
    srvx = u(2);
    srT = u(3);
    
    rho = exp(r);
    sr = sqrt(rho);
    sr1 = 1/sr;

    vx = srvx*sr1;
    T = srT*sr1;
        
    drdx  = -q(1);
    drvxdx = -q(2);
    drTdx = -q(3);
    
    dvxdx = sr1*drvxdx - 0.5*drdx*vx;
    dTdx = sr1*drTdx - 0.5*drdx*T;

    % Viscosity
    expmu = mu(12);
    expkappa = mu(13);
    nuEddy = mu(14);
    alphaEddy = mu(15);
    
    mustar = T^expmu;
%     k0 = ThermalConductivity(x,mu,eta);
    k0 = 1;
    kstar = k0*T^expkappa;
    nu = (mustar*sr1 + sr*nuEddy)/sqrt(gam*Gr);
    fc = (kstar*sr1 + sr*alphaEddy)*mw*sqrt(gam/Gr)/Pr;
    
    trr = nu*c23*2*dvxdx - 2*c23*vx/x1;
    
    fv = [0; -trr; -fc*dTdx];
end

function fv = fluxViscousWall(u,q,x,mu,eta)
    x1 = x(1);
    
    gam = mu(1);
    Gr = mu(2);
    Pr = mu(3);
    c23 = 2.0/3.0;
    mw = weightedMass(x,mu,eta);
    
    r = u(1);    
    rho = exp(r);
    sr = sqrt(rho);
    sr1 = 1/sr;
    vx = 0.0;
    T = 1.0;
        
    drdx  = -q(1);
    drvxdx = -q(2);
    drTdx = -q(3);
    
    dvxdx = sr1*drvxdx - 0.5*drdx*vx;
    dTdx = sr1*drTdx - 0.5*drdx*T;

    % Viscosity
    expmu = mu(12);
    expkappa = mu(13);
    nuEddy = mu(14);
    alphaEddy = mu(15);
    
    mustar = T^expmu;
%     k0 = ThermalConductivity(x,mu,eta);
    k0 = 1;
    kstar = k0*T^expkappa;
    nu = (mustar*sr1 + sr*nuEddy)/sqrt(gam*Gr);
    fc = (kstar*sr1 + sr*alphaEddy)*mw*sqrt(gam/Gr)/Pr;
    
    trr = nu*c23*2*dvxdx - 2*c23*vx/x1;
    
    fv = [0; -trr; -fc*dTdx];
end

function   tau = gettau(uhat, mu, eta, x, n)
    gam = mu(1);
    Gr = mu(2);
    Pr = mu(3);
    mw = weightedMass(x,mu,eta);
    
    r = uhat(1);
    srvx = uhat(2);
    srT = uhat(3);
    
    rho = exp(r);
    sr = sqrt(rho);
    sr1 = 1/sr;
    T = srT*sr1;

    vx = srvx*sr1;
    c = sqrt(T);
    
%     tauA = sqrt(vx*vx) + c;
    tauA = mu(22);
    
    % Viscosity
    expmu = mu(12);
    expkappa = mu(13);
    nuEddy = mu(14);
    alphaEddy = mu(15);
    
    mustar = T^expmu;
%     k0 = ThermalConductivity(x,mu,eta);
    k0 = 1;
    kstar = k0*T^expkappa;
    tauDv = (mustar*sr1 + sr*nuEddy)/sqrt(gam*Gr);
    tauDT = (kstar*sr1 + sr*alphaEddy)*mw*sqrt(gam/Gr)/Pr;
    
    tau(1) = tauA;
    tau(2) = tauA + tauDv;
    tau(3) = tauA + tauDT; 
end

function kappa = ThermalConductivity(x, mu, eta)
    %Parameters that need to be defined
    nspecies = 4;
    nWaves = 37;
        
    %Position
    r = x(1);
    
    R0 = mu(16);
    H0 = mu(18);
    
    %Compute weighted density compositions (n_i/rho = Chi/mi)
    Chi = sym(zeros(nspecies,1));
    Chi(1) = 1.0;
    for iSpecies = 2:nspecies
        coeffsDensity = eta((3+nspecies)*nWaves+4*(iSpecies-2)+1:(3+nspecies)*nWaves+4*(iSpecies-1));
        Chi(iSpecies) = coeffsDensity(1)*exp(coeffsDensity(2)*(r-R0)*H0) + coeffsDensity(3)*exp(coeffsDensity(4)*(r-R0)*H0);
        Chi(1) = Chi(1) - Chi(iSpecies);
    end
        
    %Compute thermal conductivity
    kappa = sym(0.0);
    ckappai = eta((3+nspecies)*nWaves+4*(nspecies-1)+nspecies+1:(3+nspecies)*nWaves+4*(nspecies-1)+2*nspecies);
    
    for iSpecies = 1:nspecies
        kappa = kappa + ckappai(iSpecies)*Chi(iSpecies);
    end
end

function [mw,dmdr] = weightedMass(x,mu,eta)
    nspecies = 4;
    nWaves = 37;

    %Position
    r = x(1);

    R0 = mu(16);
    H0 = mu(18);

    %Compute weighted density compositions (n_i/rho = Chi/mi)
    Chi = sym(zeros(nspecies,1));
    dChidr = sym(zeros(nspecies,1));
    Chi(1) = 1.0;
    for iSpecies = 2:nspecies
        coeffsDensity = eta((3+nspecies)*nWaves+4*(iSpecies-2)+1:(3+nspecies)*nWaves+4*(iSpecies-1));
        Chi(iSpecies) = coeffsDensity(1)*exp(coeffsDensity(2)*(r-R0)*H0) + coeffsDensity(3)*exp(coeffsDensity(4)*(r-R0)*H0);
        Chi(1) = Chi(1) - Chi(iSpecies);

        dChidr(iSpecies) = (coeffsDensity(1)*coeffsDensity(2)*exp(coeffsDensity(2)*(r-R0)*H0) + coeffsDensity(3)*coeffsDensity(4)*exp(coeffsDensity(4)*(r-R0)*H0))*H0;
        dChidr(1) = dChidr(1) - dChidr(iSpecies);
    end

    mass = eta((3+nspecies)*nWaves+4*(nspecies-1)+1:(3+nspecies)*nWaves+4*(nspecies-1)+nspecies);
    mw = sym(0.0);
    dmdr = sym(0.0);
    for iSpecies = 1:nspecies
        mw = mw + mass(iSpecies)*Chi(iSpecies);
        dmdr = dmdr + mass(iSpecies)*dChidr(iSpecies);
    end
end