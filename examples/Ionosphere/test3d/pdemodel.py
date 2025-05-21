from numpy import array
from sympy import sqrt

# mu = [m, kB, e, g]
# v = [E1, E2, E3, B1, B2, B3, vn1, vn2, vn3, nu, P, L, H, C]

def mass(u, q, w, v, x, t, mu, eta):
    m = array([1, 1, 1, 1, 1])
    return m

def flux(u, q, w, v, x, t, mu, eta):
    m = mu[0]
    kB = mu[1]
    u1 = u[0]
    u2 = u[1:4]
    u3 = u[4]
    kappa = 1.24e4*(u3/(kB*u1))**2.5
    q1 = array([q[0], q[5], q[10]])
    q3 = array([q[4], q[9], q[14]])
    f1 = u2/m
    f21 = [u2[0]**2/(m*u1) + u3, u2[0]*u2[1]/(m*u1)  , u2[0]*u2[2]/(m*u1)  ]
    f22 = [u2[1]*u2[0]/(m*u1)  , u2[1]**2/(m*u1) + u3, u2[1]*u2[2]/(m*u1)  ]
    f23 = [u2[2]*u2[0]/(m*u1)  , u2[2]*u2[1]/(m*u1)  , u2[2]**2/(m*u1) + u3]
    f3 = u3*u2/(m*u1) - 2*kappa*(u3*q1 - u1*q3)/(3*kB*u1**2)
    f = array([f1[0], f21[0], f22[0], f23[0], f3[0],
               f1[1], f21[1], f22[1], f23[1], f3[1],
               f1[2], f21[2], f22[2], f23[2], f3[2]])
    return f

def source(u, q, w, v, x, t, mu, eta):
    m = mu[0]
    e = mu[2]
    g = mu[3]
    E = v[0:3]
    B = v[3:6]
    vn = v[6:9]
    nu = v[9]
    P = v[10]
    L = v[11]
    H = v[12]
    C = v[13]
    u1 = u[0]
    u2 = u[1:4]
    u3 = u[4]
    q1 = [q[0], q[5], q[10]]
    trq2 = q[1] + q[7] + q[13]
    r = sqrt(x[0]**2 + x[1]**2 + x[2]**2)

# start with a simple case
    nu = 0
    P = 0
    L = 0
    H = 0
    C = 0

    s1 = P - L*u1
    s21 = (P/u1-L)*u2[0] + u1*m*g*x[0]/r + u1*e*E[0] + e*(u2[1]*B[2]-u2[2]*B[1])/m + nu*(u1*m*vn[0]-u2[0])
    s22 = (P/u1-L)*u2[1] + u1*m*g*x[1]/r + u1*e*E[1] + e*(u2[2]*B[0]-u2[0]*B[2])/m + nu*(u1*m*vn[1]-u2[1])
    s23 = (P/u1-L)*u2[2] + u1*m*g*x[2]/r + u1*e*E[2] + e*(u2[0]*B[1]-u2[1]*B[0])/m + nu*(u1*m*vn[2]-u2[2])
    s3 = 2*u3*(u1*trq2 - (u2[0]*q1[0]+u2[1]*q1[1]+u2[2]*q1[2]))/(3*m*u1) + u3*(P/u1 - L) + 2*(H - C*u3)/3

    s = array([s1, s21, s22, s23, s3])
    return s

def ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    m = mu[0]
    vn = v[6:9]
    P = v[10]
    L = v[11]
    H = v[12]
    C = v[13]
    ub = array([P/L, m*P*vn[0]/L, m*P*vn[1]/L, m*P*vn[2]/L, H/C])
    return [ub, u]

def fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    f = flux(u, q, w, v, x, t, mu, eta)
    fb = array([f[:, 0]*n[0] + f[:, 1]*n[1] + f[:, 2]*n[2] + tau*(u-uhat)])

    fbu = array([0, 0, 0, 0, 0])
    return [fb, fbu]

def fbouhdg(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    f = flux(u, q, w, v, x, t, mu, eta)
    fb = array([f[:, 0]*n[0] + f[:, 1]*n[1] + f[:, 2]*n[2] + tau*(u-uhat)])

    m = mu[0]
    vn = v[6:9]
    P = v[10]
    L = v[11]
    H = v[12]
    C = v[13]
    ub = array([P/L, m*P*vn[0]/L, m*P*vn[1]/L, m*P*vn[2]/L, H/C])
    fl = tau*(ub - uhat)

    return [fl, fb]

def initu(x, mu, eta):
    u0 = array([0, 0, 0, 0, 0])
    return u0
