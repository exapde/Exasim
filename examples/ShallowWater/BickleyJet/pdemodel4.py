from numpy import array, reshape
from sympy import sin, pi, exp, sqrt, tan, cos, cosh

def mass(u, q, w, v, x, t, mu, eta):
    m = array([1.0, 1.0]);
    return m;

def flux(u, q, w, v, x, t, mu, eta):
    r = w[0];  # pressure
    uv = u[0];
    vv = u[1];    
    p = r;        
    f = array([uv*uv+p, vv*uv, uv*vv, vv*vv+p]);
    f = reshape(f,(2,2),'F');    
    return f;

def source(u, q, w, v, x, t, mu, eta):
    s = array([0.0, 0.0]);
    return s;

def sourcew(u, q, w, v, x, t, mu, eta):
    s = array([mu[0]*(q[0]+q[3])]);
    return s;

def fbou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    f = flux(u, q, w, v, x, t, mu, eta);
    fb = f[:,0]*n[0] + f[:,1]*n[1] + tau[0]*(u-uhat);
    return fb;

def ubou(u, q, w, v, x, t, mu, eta, uhat, n, tau):
    ub = array([0.0, 0.0]);
    return ub;

def initu(x, mu, eta):
    epsil = 0.1; # perturbation magnitude
    l = 0.5;     # Gaussian width
    k = 0.5;     # Sinusoidal wavenumber
    
    x1 = x[0];
    x2 = x[1];
    
    # The Bickley jet
    U = (1/cosh(x2))*(1/cosh(x2));

    # Slightly off-center vortical perturbations
    Psiprime = exp(-(x2 + l/10)*(x2 + l/10) / (2*(l*l))) * cos(k * x1) * cos(k * x2);

    # Vortical velocity fields (ũ, ṽ) = (-∂_y, +∂_x) ψ̃
    uprime =  Psiprime * (k * tan(k * x2) + x2 /(l*l)); 
    vprime = -Psiprime * k * tan(k * x1); 

    u01 = U + epsil * uprime; # h*u
    u02 = epsil * vprime;  # h*v
   
    u0 = array([u01, u02]);
    return u0;

def initw(x, mu, eta):

    w0 = array([0.0]);
    return w0;

