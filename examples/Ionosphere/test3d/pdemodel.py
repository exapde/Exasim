from numpy import array, reshape, hstack, zeros
from sympy import sqrt, exp, log

def mass(u, q, w, v, x, t, mu, eta):
    m = array([1.0, 1.0, 1.0, 1.0, 1.0])
    return m

def flux(u, q, w, v, x, t, mu, eta):
    gam = mu[0]
    
    ln = u[0]
    u1 = u[1]
    u2 = u[2]
    u3 = u[3]
    T  = u[4]
    p = T/gam

    lnx = -q[0]
    u1x = -q[1]
    u2x = -q[2]
    u3x = -q[3]
    Tx  = -q[4]

    lny = -q[5]
    u1y = -q[6]
    u2y = -q[7]
    u3y = -q[8]
    Ty  = -q[9]

    lnz = -q[10]
    u1z = -q[11]
    u2z = -q[12]
    u3z = -q[13]
    Tz  = -q[14]

    fi = array([ln*u1, u1*u1+p, u1*u2, u1*u3, T*u1,
                ln*u2, u1*u2, u2*u2+p, u2*u3, T*u2,
                ln*u3, u1*u3, u2*u3, u3*u3+p, T*u3])
    
    f = fi
    f = reshape(f,(5,3),'F');
    return f

def source(u, q, w, v, x, t, mu, eta):
    gam = mu[0]
    Ro = mu[1]

    R0 = mu[8]
    x1 = x[0]
    x2 = x[1]
    x3 = x[2]
    r  = sqrt(x1**2 + x2**2 + x3**2)
    g0 = 1/gam
    g = g0*(R0/r)**2

    ln = u[0]
    u1 = u[1]
    u2 = u[2]
    u3 = u[3]
    T  = u[4]
    p = T/gam

    lnx = -q[0]
    u1x = -q[1]
    u2x = -q[2]
    u3x = -q[3]
    Tx  = -q[4]

    lny = -q[5]
    u1y = -q[6]
    u2y = -q[7]
    u3y = -q[8]
    Ty  = -q[9]

    lnz = -q[10]
    u1z = -q[11]
    u2z = -q[12]
    u3z = -q[13]
    Tz  = -q[14]

    # acceleration
    ax = -g*x1/r + x1/Ro**2 + 2*u2/Ro
    ay = -g*x2/r + x2/Ro**2 - 2*u1/Ro
    az = -g*x3/r

    # divergence of velocity
    divV = u1x + u2y + u3z

    # Lorentz force: E + vxB
    piE = mu[2]
    piB = mu[3]
    Bx = v[0]
    By = v[1]
    Bz = v[2]
    Ex = v[3]
    Ey = v[4]
    Ez = v[5]
    fLx = piE*Ex + piB*(u2*Bz - u3*By)
    fLy = piE*Ey + piB*(u3*Bx - u1*Bz)
    fLz = piE*Ez + piB*(u1*By - u2*Bx)

    # source terms
    sln = (ln-1)*divV
    su1 = ax + divV*u1 - p*lnx + fLx
    su2 = ay + divV*u2 - p*lny + fLy
    su3 = az + divV*u3 - p*lnz + fLz
    sT  = (2 - gam)*T*divV

    s = array([sln, su1, su2, su3, sT])
    return s

def ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    ub = zeros((5,2));
    return ub

def fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    fb = zeros((5,2));
    return fb

def fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    fl = 0*u
    fl[0] = 0.0 - uhat[0]
    fl[1] = 0.0 - uhat[1]
    fl[2] = 0.0 - uhat[2]
    fl[3] = 0.0 - uhat[3]
    fl[4] = 1.0 - uhat[4]

    fu = 0*u
    fu[0] = u[0] - uhat[0]
    fu[1] = 0.0 - uhat[1]
    fu[2] = 0.0 - uhat[2]
    fu[3] = 0.0 - uhat[3]
    fu[4] = 5.0 - uhat[4]

    fb = hstack((fl, fu)); # wall and freestream boundary conditions
    fb = reshape(fb,(5,2),'F');
    return fb

def initu(x, mu, eta):
    x1 = x[0]
    x2 = x[1]
    x3 = x[2]
    r = sqrt(x1**2 + x2**2 + x3**2)

    Tbot = 1.0
    Ttop = 5.0
    R0 = mu[8]
    H0 = mu[10]
    h0 = 40000.0/H0
    a0 = -1

    T = Ttop - (Ttop-Tbot)*exp(-(r-R0)/h0)
    logp_p0 = a0*h0/Ttop*log(1+Ttop/Tbot*(exp((r-R0)/h0)-1))
    ln = logp_p0 - log(T)


    u0 = array([ln, 0.0, 0.0, 0.0, T])
    return u0
