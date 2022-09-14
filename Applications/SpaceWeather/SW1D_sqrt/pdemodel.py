"""
Module that includes GITM 1D (reduced stiffness formulation) model functions.

Authors: Opal Issan
Version: Sept 9, 2022
"""
from numpy import array, arange
from copy import deepcopy
from sympy import exp, sqrt, log, pi, sin , cos, tanh


def mass(u, q, w, v, x, t, mu, eta):
    m = array([1.0, 1.0, 1.0])
    return m


def flux(u, q, w, v, x, t, mu, eta):
    fi = fluxInviscid(u, mu)
    fv = fluxViscous(u, q, x, mu)
    return fi + fv


def fluxWall(u, q, w, v, x, t, mu, eta):
    fi = fluxInviscidWall(u, mu)
    fv = fluxViscousWall(u, q, x, mu)
    return fi + fv


def source(u, q, w, v, x, t, mu, eta):
    x1 = x[0]
    c0 = mu[13]
    gam = mu[0]
    gam1 = gam - 1
    Gr = mu[1]
    Pr = mu[2]
    c23 = 2.0 / 3.0

    r = u[0]
    srvx = u[1]
    srT = u[2]

    rho = exp(r)
    sr = sqrt(rho)
    sr1 = 1 / sr
    r_1 = r - 1

    vx = srvx * sr1
    T = srT * sr1
    p = srT / gam

    drdx = -q[0]
    drvxdx = -q[1]
    drTdx = -q[2]

    dvxdx = sr1 * drvxdx - 0.5 * drdx * vx
    dTdx = sr1 * drTdx - 0.5 * drdx * T

    # Viscosity
    mustar = sqrt(T)
    kstar = T ** 0.75
    nu = mustar * sr1 / sqrt(gam * Gr)
    fc = kstar * sr1 * sqrt(gam / Gr) / Pr

    trr = nu * c23 * 2 * dvxdx - c23 * c0 * vx / x1
    trd = nu * 2 * c0 * (dvxdx - vx / x1) / x1

    R0 = mu[9]
    gravity0 = 1 / gam
    gravity = gravity0 * R0 ** 2 / (x1 ** 2)
    Fr = mu[3]
    ar = -gravity + x1 * Fr ** 2 / gam

    trp = 2 * c23 * nu * (dvxdx ** 2 - c0 * vx * dvxdx / x1 + 0.5 * c0 * (3 - c0) * vx ** 2 / x1 ** 2)
    SigmadV = gam * gam1 * trp

    q_EUV = EUVsource1D(u, x, t, mu, eta)

    return array([r_1 * dvxdx - c0 * vx / x1,
                  sr * ar + 0.5 * (dvxdx - c0 * vx / x1) * srvx - 0.5 * p * drdx + 0.5 * trr * drdx + 0.5 * trd,
                  sr * q_EUV + (3 / 2 - gam) * srT * dvxdx + c0 * (1 / 2 - gam) * srT * vx / x1 +
                  fc * dTdx * (c0 / x1 + 0.5 * drdx) + SigmadV])


def fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    tau = gettau(uhat, mu, eta, n)

    f = fluxWall(u, q, w, v, x, t, mu, eta)
    fw = f * n  # numerical flux at freestream boundary
    fw[0] = fw[0] + tau[0] * (u[0] - uhat[0])
    fw[1] = fw[1] + tau[1] * (u[1] - uhat[1])
    fw[2] = fw[2] + tau[2] * (u[2] - uhat[2])

    # Inviscid outer boundary
    fi2 = fluxInviscid(u, mu)
    ft = fi2 * n
    ft[0] = ft[0] + tau[0] * (u[0] - uhat[0])
    ft[1] = ft[1] + tau[1] * (u[1] - uhat[1])
    ft[2] = ft[2] + tau[2] * (u[2] - uhat[2])

    return array([fw, ft])


def ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    Tbot = mu[7]

    # Isothermal Wall
    r = u[0]
    rho = exp(r)
    sr = sqrt(rho)

    utw1 = u
    utw1[1] = 0.0
    utw1[2] = sr * Tbot

    # Inviscid wall
    utw2 = u

    return array([utw1, utw2])


def initu(x, mu, eta):
    x1 = x[0]

    Fr = mu[3]

    Tbot = mu[7]
    Ttop = mu[8]
    R0 = mu[9]
    Ldim = mu[11]
    h0 = 35000 / Ldim

    a0 = (-1 + Fr ** 2 * R0)

    T = Ttop - (Ttop - Tbot) * exp(-(x1 - R0) / h0)
    logp_p0 = a0 * h0 / Ttop * log(1 + Ttop / Tbot * (exp((x1 - R0) / h0) - 1))
    rtilde = logp_p0 - log(T)
    rho = exp(rtilde)
    srT = sqrt(rho) * T

    return array([rtilde, 0.0, srT])


def stab(u1, q1, w1, v1, x, t, mu, eta, uhat, n, tau, u2, q2, w2, v2):
    uhat = 0.5 * (u1 + u2)
    tau = gettau(uhat, mu, eta, n)
    ftau = tau
    ftau[0] = tau[0] * (u1[0] - u2[0])
    ftau[1] = tau[1] * (u1[1] - u2[1])
    ftau[2] = tau[2] * (u1[2] - u2[2])
    return ftau


def fluxInviscid(u, mu):
    gam = mu[0]
    r = u[0]
    srvx = u[1]
    srT = u[2]

    rho = exp(r)
    sr = sqrt(rho)
    sr1 = 1 / sr

    vx = srvx * sr1

    p = srT / gam

    return array([r * vx,
                  srvx * vx + p,
                  srT * vx])


def fluxInviscidWall(u, mu):
    gam = mu[0]

    r = u[0]
    rho = exp(r)
    sr = sqrt(rho)

    Tbot = mu[7]
    p = sr * Tbot / gam

    return array([0, p, 0])


def fluxViscous(u, q, x, mu):
    x1 = x[0]
    c0 = mu[13]

    gam = mu[0]
    Gr = mu[1]
    Pr = mu[2]
    c23 = 2.0 / 3.0

    r = u[0]
    srvx = u[1]
    srT = u[2]

    rho = exp(r)
    sr = sqrt(rho)
    sr1 = 1 / sr

    vx = srvx * sr1
    T = srT * sr1

    drdx = -q[0]
    drvxdx = -q[1]
    drTdx = -q[2]

    dvxdx = sr1 * drvxdx - 0.5 * drdx * vx
    dTdx = sr1 * drTdx - 0.5 * drdx * T

    # Viscosity
    mustar = sqrt(T)
    kstar = T ** 0.75
    nu = mustar * sr1 / sqrt(gam * Gr)
    fc = kstar * sr1 * sqrt(gam / Gr) / Pr

    trr = nu * c23 * 2 * dvxdx - c23 * c0 * vx / x1

    return array([0, -trr, -fc * dTdx])


def fluxViscousWall(u, q, x, mu):
    x1 = x[0]
    c0 = mu[13]

    gam = mu[0]
    Gr = mu[1]
    Pr = mu[2]
    c23 = 2.0 / 3.0

    r = u[0]
    rho = exp(r)
    sr = sqrt(rho)
    sr1 = 1 / sr
    vx = 0.0
    T = mu[7]

    drdx = -q[0]
    drvxdx = -q[1]
    drTdx = -q[2]

    dvxdx = sr1 * drvxdx - 0.5 * drdx * vx
    dTdx = sr1 * drTdx - 0.5 * drdx * T

    # Viscosity
    mustar = sqrt(T)
    kstar = T ** 0.75
    nu = mustar * sr1 / sqrt(gam * Gr)
    fc = kstar * sr1 * sqrt(gam / Gr) / Pr

    trr = nu * c23 * 2 * dvxdx - c23 * c0 * vx / x1

    return array([0, -trr, -fc * dTdx])


def gettau(uhat, mu, eta, n):
    gam = mu[0]
    Gr = mu[1]
    Pr = mu[2]

    r = uhat[0]
    # srvx = uhat[1]
    srT = uhat[2]

    rho = exp(r)
    sr = sqrt(rho)
    sr1 = 1 / sr
    T = srT * sr1

    # vx = srvx * sr1
    # c = sqrt(T)
    #  tauA = sqrt(vx * vx) + c
    tauA = mu[17]

    # Viscosity
    mustar = sqrt(T)
    kstar = T ** 0.75
    tauDv = mustar * sr1 / sqrt(gam * Gr)
    tauDT = kstar * sr1 * sqrt(gam / Gr) / Pr

    return array([tauA, tauA + tauDv, tauA + tauDT])


def EUVsource1D(u, x, t, mu, eta):
    r = x[0]

    gam = mu[0]
    gam1 = gam - 1.0

    Fr = mu[3]
    omega = Fr / sqrt(gam)
    K0 = mu[4]
    M0 = mu[5]

    R0 = mu[9]

    longitude = mu[14] * pi / 180
    latitude = mu[15] * pi / 180
    declinationSun = mu[16] * pi / 180

    # computation of angles
    # define local time
    long_offset = omega * t - pi / 2
    localTime = longitude + long_offset
    cosChi = sin(declinationSun) * sin(latitude) + cos(declinationSun) * cos(latitude) * cos(localTime)
    # cosChi = cos(localTime)

    absSinChi = sqrt(1 - cosChi ** 2)

    # Computation F10.7(let's assume it constant at first, the variation is at another scale)
    F10p7 = 100
    F10p7_81 = 100
    F10p7_mean = 0.5 * (F10p7 + F10p7_81)

    rtilde = u[0]
    rho = exp(rtilde)
    T = u[2] / sqrt(rho)

    # Quantities
    gravity = (R0 ** 2 / (r ** 2)) / gam
    H = T / (gam * gravity)

    # Chapman integral
    Rp = rho * H
    Xp = r / H
    y = sqrt(Xp / 2) * abs(cosChi)

    Ierf = 0.5 * (1 + tanh(1000 * (8 - y)))
    a_erf = 1.06069630
    b_erf = 0.55643831
    c_erf = 1.06198960
    d_erf = 1.72456090
    f_erf = 0.56498823
    g_erf = 0.06651874

    erfcy = Ierf * (a_erf + b_erf * y) / (c_erf + d_erf * y + y * y) + (1 - Ierf) * f_erf / (g_erf + y)

    IcosChi = 0.5 * (1 + tanh(100000 * cosChi))
    IsinChi = 0.5 * (1 + tanh(100000 * (r * absSinChi - R0)))

    alpha1 = Rp * erfcy * sqrt(0.5 * pi * Xp)
    auxXp = (1 - IcosChi) * IsinChi * Xp * (1 - absSinChi)
    Rg = rho * H * exp(auxXp)
    alpha2 = (2 * Rg - Rp * erfcy) * sqrt(0.5 * pi * Xp)

    alpha = IcosChi * alpha1 + (1 - IcosChi) * (IsinChi * alpha2 + (1 - IsinChi) * 1e32)

    Q = 0
    for iWave in arange(0, 37):
        lambda_1 = eta[iWave]
        crossSection = eta[37+iWave]
        AFAC = eta[2 * 37+iWave]
        F74113 = eta[3 * 37+iWave]
        tau = M0 * crossSection * alpha

        slope0 = 1 + AFAC * (F10p7_mean-80)
        Islope = 0.5 * (1+tanh(1000 * (slope0-0.8)))
        slopeIntensity = slope0 * Islope + 0.8 * (1-Islope)
        Intensity0 = F74113 * slopeIntensity
        Intensity = Intensity0 * exp(-tau)

        Q = Q + crossSection * Intensity / lambda_1

    eff = mu[12]
    return gam * gam1 * eff * Q / K0