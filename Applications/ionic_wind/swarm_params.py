import numpy as np

def get_swarm_params(E, N):
    # Using the values in the appendix of this paper: https://iopscience.iop.org/article/10.1088/0022-3727/30/4/017/pdf
    # Note that the reduced E field is converted to V*cm2 for performing the calculations, then the resultant quantity is converted back from cm->m.

    cm2m = 0.01
    m2cm = 100
    E_N = E/N*m2cm**2   # Reduced electric field in V*cm2. Expects E to be the L2 norm of the E field.

    # First ionization coefficient
    if E_N <= 1.5e-15:
        alpha = 6.619e-17*np.exp(-5.593e-15/E_N)*cm2m**2

    elif E_N > 1.5e-15:
        alpha = 2e-16*np.exp(-7.248e-15/E_N)*cm2m**2

    # Two-body attachment coefficient
    if E_N <= 1.05e-15:
        eta2 = 6.089e-4*E_N - 2.893e-19*cm2m**2

    elif E_N > 1.05e-15:
        eta2 = 8.889e-5*E_N + 2.567e-19*cm2m**2

    # Recombination coefficient
    beta = 2e-7*cm2m**3

    # Electron mobility
    # Removing the sign term because these formulas were designed for computing the velocity, not the mobility. The mobility always has the same sign.

    if E_N <= 2.6e-17:
        mue = -(3.38e4 + 6.87e22*E_N) * cm2m / E   # Formula is for the velocity, so need to divide by the E field strength to get the mobility
    
    elif E_N > 2.6e-17 and E_N <= 1e-16:
        mue = -(1.63e6 + 7.2973e21*E_N) * cm2m / E   # Formula is for the velocity, so need to divide by the E field strength to get the mobility

    elif E_N > 1e-16 and E_N <= 2e-15:
        mue = -(1.3e6 + 1.03e22*E_N) * cm2m / E   # Formula is for the velocity, so need to divide by the E field strength to get the mobility

    elif E_N > 2e-15:
        mue = -(7.1e6 + 7.4e21*E_N) * cm2m / E   # Formula is for the velocity, so need to divide by the E field strength to get the mobility
    
    # Negative ion mobility
    if E_N <= 5e-16:
        mun = -1.86*cm2m**2
    elif E_N > 5e-16:
        mun = -2.7*cm2m**2

    # Positive ion mobility
    # Assuming that P0/P is approximately 1 and that the units for the E field are V/cm
    mup = 2.34*cm2m**2

    # Electron diffusion coefficient
    # Unit conversion (cm->m) is "baked in" through the mue calculated previously
    D = mue*0.3341e9*E_N**0.54069

    return alpha, eta2, beta, D, mue, mup, mun


"""
Calc the neutral number density from the ideal gas law
PV=Nk_bT    https://en.wikipedia.org/wiki/Ideal_gas_law


"""

P = 101325 # Pa
V = 1 # m^3
T = 273.15 # K
k_b = 1.380649e-23 #m2 kg s-2 K-1

N = P*V/(k_b*T)
E_bd = 3e6
print(N)
alpha, eta2, beta, D, mue, mup, mun = get_swarm_params(3e6, N)

print(f'mue at E_bd: {mue}')
print(f'D at E_bd: {D}')