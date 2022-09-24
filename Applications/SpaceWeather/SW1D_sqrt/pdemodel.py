from numpy import array, reshape, hstack
from sympy import exp, sqrt, log, pi, tanh, sin, cos

def mass(u, q, w, v, x, t, mu, eta):
    m = array([1.0, 1.0, 1.0])
    return m

def flux(u, q, w, v, x, t, mu, eta):
    fi = fluxInviscid(u,mu)
    fv = fluxViscous(u,q,x,mu)
      
    f = fi+fv
    return f

def fluxWall(u, q, w, v, x, t, mu, eta):
    fi = fluxInviscidWall(u,mu)
    fv = fluxViscousWall(u,q,x,mu)
      
    f = fi+fv
    return f

def source(u, q, w, v, x, t, mu, eta):
    x1 = x[0]
    c0 = mu[13]

    gam = mu[0]
    gam1 = gam-1
    Gr = mu[1]
    Pr = mu[2]
    c23 = 2.0/3.0
    
    r = u[0]
    srvx = u[1]
    srT = u[2]
    
    rho = exp(r)
    sr = sqrt(rho)
    sr1 = 1/sr
    r_1 = r-1
    
    vx = srvx*sr1
    T = srT*sr1
    p = srT/gam
    
    drdx  = -q[0]
    drvxdx = -q[1]
    drTdx = -q[2]
        
    dvxdx = sr1*drvxdx - 0.5*drdx*vx
    dTdx = sr1*drTdx - 0.5*drdx*T
    
    # Viscosity
    mustar = sqrt(T)
    kstar = T**0.75
    nu = mustar*sr1/sqrt(gam*Gr)
    fc = kstar*sr1*sqrt(gam/Gr)/Pr
    
    trr = nu*c23*2*dvxdx - c23*c0*vx/x1
    trd = nu*2*c0*(dvxdx-vx/x1)/x1
    
    R0 = mu[9]
    gravity0 = 1/gam
    gravity = gravity0*R0**2/x1**2
    Fr = mu[3]
    ar = -gravity + x1*Fr**2/gam
    
    trp = 2*c23*nu*(dvxdx**2 - c0*vx*dvxdx/x1 + 0.5*c0*(3-c0)*vx**2/(x1**2))
    SigmadV = gam*gam1*trp
    
    q_EUV = EUVsource1D(u, x, t, mu, eta)
    
    s = array([r_1*dvxdx - c0*vx/x1,
        sr*ar + 0.5*(dvxdx-c0*vx/x1)*srvx - 0.5*p*drdx + 0.5*trr*drdx + 0.5*trd,
        sr*q_EUV + (3/2-gam)*srT*dvxdx + c0*(1/2-gam)*srT*vx/x1 + fc*dTdx*(c0/x1 + 0.5*drdx) + SigmadV])
        
    return s

def fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    #tau = gettau(uhat, mu, eta, n)
    tau = array([0.0, 0.0, 0.0])

    f = fluxWall(u, q, w, v, x, t, mu, eta)
    fw0 = f[0]*n[0] + tau[0]*(u[0]-uhat[0])
    fw1 = f[1]*n[0] + tau[1]*(u[1]-uhat[1])
    fw2 = f[2]*n[0] + tau[2]*(u[2]-uhat[2])
    fw = array([fw0, fw1, fw2])
    
    # Inviscid outer boundary
    fi = fluxInviscid(u,mu)
    ft0 = fi[0]*n[0] + tau[0]*(u[0]-uhat[0])
    ft1 = fi[1]*n[0] + tau[1]*(u[1]-uhat[1])
    ft2 = fi[2]*n[0] + tau[2]*(u[2]-uhat[2])
    ft = array([ft0, ft1, ft2])

    fb = hstack((fw, ft))
    fb = reshape(fb,[3,2],'F')
    return fb

def ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    Tbot = mu[7]

    # Isothermal Wall
    r = u[0]
    rho = exp(r)
    sr = sqrt(rho)
    
    utw1 = array([r, 0.0, sr*Tbot])

    # Inviscid wall
    utw2 = u
        
    ub = hstack((utw1, utw2)) 
    ub = reshape(ub,(3,2),'F')
    return ub

def initu(x, mu, eta):
    x1 = x[0]
    
    Fr = mu[3]

    Tbot = mu[7]
    Ttop = mu[8]
    R0 = mu[9]
    Ldim = mu[11]
    h0 = 35000/Ldim

    a0 = (-1 + (Fr**2)*R0)
    
    T = Ttop - (Ttop-Tbot)*exp(-(x1-R0)/h0)
    logp_p0 = a0*h0/Ttop*log(1+Ttop/Tbot*(exp((x1-R0)/h0)-1))
    rtilde = logp_p0 - log(T)
    rho = exp(rtilde)
    srT = sqrt(rho)*T
    
    u0 = array([rtilde, 0.0, srT])
    return u0


def stab(u1, q1, w1, v1, x, t, mu, eta, uhat, n, tau, u2, q2, w2, v2):
    uhat = 0.5*(u1+u2)
    tau = gettau(uhat, mu, eta, n)
    
    ftau0 = tau[0]*(u1[0] - u2[0])
    ftau1 = tau[1]*(u1[1] - u2[1])
    ftau2 = tau[2]*(u1[2] - u2[2])
    ftau = array([ftau0, ftau1, ftau2])
    return ftau


def fluxInviscid(u,mu):
    gam = mu[0]  
    r = u[0]
    srvx = u[1]
    srT = u[2]
    
    rho = exp(r)
    sr = sqrt(rho)
    sr1 = 1/sr
    
    vx = srvx*sr1
    
    p = srT/gam

    fi = array([r*vx, srvx*vx+p, srT*vx])
    return fi

def fluxInviscidWall(u,mu):
    gam = mu[0]  
    r = u[0]
    rho = exp(r)
    sr = sqrt(rho)
    
    Tbot = mu[7]
    p = sr*Tbot/gam

    fi = array([0.0, p, 0.0])
    return fi

def fluxViscous(u,q,x,mu):
    x1 = x[0]
    c0 = mu[13]
    
    gam = mu[0]
    Gr = mu[1]
    Pr = mu[2]
    c23 = 2.0/3.0
    
    r = u[0]
    srvx = u[1]
    srT = u[2]
    
    rho = exp(r)
    sr = sqrt(rho)
    sr1 = 1/sr

    vx = srvx*sr1
    T = srT*sr1
        
    drdx  = -q[0]
    drvxdx = -q[1]
    drTdx = -q[2]
    
    dvxdx = sr1*drvxdx - 0.5*drdx*vx
    dTdx = sr1*drTdx - 0.5*drdx*T

    # Viscosity
    mustar = sqrt(T)
    kstar = T**0.75
    nu = mustar*sr1/sqrt(gam*Gr)
    fc = kstar*sr1*sqrt(gam/Gr)/Pr
    
    trr = nu*c23*2*dvxdx - c23*c0*vx/x1
    
    fv = array([0, -trr, -fc*dTdx])
    return fv

def fluxViscousWall(u,q,x,mu):
    x1 = x[0]
    c0 = mu[13]
    
    gam = mu[0]
    Gr = mu[1]
    Pr = mu[2]
    c23 = 2.0/3.0
    
    r = u[0]
    rho = exp(r)
    sr = sqrt(rho)
    sr1 = 1/sr
    vx = 0.0
    T = mu[7]
        
    drdx  = -q[0]
    drvxdx = -q[1]
    drTdx = -q[2]
    
    dvxdx = sr1*drvxdx - 0.5*drdx*vx
    dTdx = sr1*drTdx - 0.5*drdx*T

    # Viscosity
    mustar = sqrt(T)
    kstar = T**0.75
    nu = mustar*sr1/sqrt(gam*Gr)
    fc = kstar*sr1*sqrt(gam/Gr)/Pr
    
    trr = nu*c23*2*dvxdx - c23*c0*vx/x1
    
    fv = array([0, -trr, -fc*dTdx])
    return fv

def gettau(uhat, mu, eta, n):
    gam = mu[0]
    Gr = mu[1]
    Pr = mu[2]
    
    r = uhat[0]
    srvx = uhat[1]
    srT = uhat[2]
    
    rho = exp(r)
    sr = sqrt(rho)
    sr1 = 1/sr
    T = srT*sr1

#     vx = srvx*sr1
#     c = sqrt(T);   
#     tauA = sqrt(vx*vx) + c
    tauA = mu[17]

    # Viscosity
    mustar = sqrt(T)
    kstar = T**0.75
    tauDv = mustar*sr1/sqrt(gam*Gr)
    tauDT = kstar*sr1*sqrt(gam/Gr)/Pr
    
    tau = array([tauA, tauA + tauDv, tauA + tauDT])
    return tau

def EUVsource1D(u, x, t, mu, eta):

    r = x[0]
    
    gam = mu[0]
    gam1 = gam - 1.0
    
    Fr = mu[3]
    omega = Fr/sqrt(gam)
    K0 = mu[4]
    M0 = mu[5]
    
    R0 = mu[9]
    
    longitude = mu[14]*pi/180
    latitude = mu[15]*pi/180
    declinationSun = mu[16]*pi/180
    
    ## computation of angles
    #define local time
    long_offset = omega*t - pi/2
    localTime = longitude + long_offset
    cosChi = sin(declinationSun)*sin(latitude) + cos(declinationSun)*cos(latitude)*cos(localTime)
    
    absSinChi = sqrt(1-cosChi**2)
    
    #Computation F10.7 (let's assume it constant at first, the variation is at another scale)
    F10p7 = 100
    F10p7_81 = 100
    F10p7_mean = 0.5*(F10p7 + F10p7_81)
    
    rtilde = u[0]
    rho = exp(rtilde)
    T = u[2]/sqrt(rho)

    # Quantities
    gravity = (R0**2/(r**2))/gam
    H = T/(gam*gravity)
    
    #Chapman integral 
    Rp = rho*H
    Xp = r/H
    y = sqrt(Xp/2)*abs(cosChi)
    
    Ierf = 0.5*(1+tanh(1000*(8-y)))
    a_erf = 1.06069630
    b_erf = 0.55643831
    c_erf = 1.06198960
    d_erf = 1.72456090
    f_erf = 0.56498823
    g_erf = 0.06651874

    erfcy = Ierf*(a_erf + b_erf*y)/(c_erf + d_erf*y + y*y) + (1-Ierf)*f_erf/(g_erf + y)
    
    IcosChi = 0.5*(1 + tanh(100*cosChi))
    IsinChi = 0.5*(1 + tanh(100*(r*absSinChi - R0)))
    
    alpha1 = Rp*erfcy*sqrt(0.5*pi*Xp)
    auxXp = (1-IcosChi)*IsinChi*Xp*(1-absSinChi)
    Rg = rho*H*exp(auxXp)
    alpha2 = (2*Rg - Rp*erfcy)*sqrt(0.5*pi*Xp)
    
    alpha = IcosChi*alpha1 + (1-IcosChi)*(IsinChi*alpha2 + (1-IsinChi)*1e2)
    
    Q = 0
    for iWave in range(0,37):
        lambdaw = eta[iWave]
        crossSection = eta[37+iWave]
        AFAC = eta[2*37+iWave]
        F74113 = eta[3*37+iWave]
        
        tau = M0*crossSection*alpha
        
        slope0 = 1 + AFAC*(F10p7_mean-80)
        Islope = 0.5*(1+tanh(1000*(slope0-0.8)))
        slopeIntensity =  slope0*Islope + 0.8*(1-Islope)
        Intensity0 = F74113*slopeIntensity
        Intensity = Intensity0*exp(-tau)
        
        Q = Q + crossSection*Intensity/lambdaw
    
    eff = mu[12]
    s_EUV = gam*gam1*eff*Q/K0

    return s_EUV