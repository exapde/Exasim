import sys
import os
mppPyDir = os.environ.get('MPP_LOCALPY')
sys.path.append(mppPyDir)
import _mutationpp as mpp
import numpy as np


def print_equil_vals(mix, T, P):
    mix.equilibrate(T, P)
    rhovec = mix.densities()
    rho = np.sum(rhovec)
    e = mix.mixtureEnergyMass()
    rhoE = rho * e + 0.0
    w = mix.netProductionRates()

    print("To check")
    print("(T, P) = (" + str(mix.T()) + ", " + str(mix.P()) + ")")
    print("Initial conditions (dimensional): ")
    print("(rho_i, rhou, rhoE) = (" + str(rhovec) + " 0.0 " + str(rhoE) + ")")
    print("omega: " + str(w))

def init_mixture():
    opt = mpp.MixtureOptions("air_5")
    opt.setStateModel("ChemNonEq1T")
    opt.setViscosityAlgorithm("Gupta-Yos")
    opt.setThermodynamicDatabase("NASA-9")

    mix = mpp.Mixture(opt)
    return mix

def rhoe(mix, u):
    e = mix.mixtureEnergyMass()
    rho = mix.density()
    uTu2 = 0.5 * u**2 # Assumes scalar density
    rhoe = rho * e - rho * uTu2
    return rhoe

# def return_farfield():
u_inf = 5956
rho_inf = 1.547e-3
# p_inf = 476
T_inf = 901
# Mach_inf = 8.98
Y_N = 6.5e-7
Y_O = 0.22830098
Y_NO = 0.01026010
Y_N2 = 0.75430704
Y_O2 = 0.00713123
Y_vec = np.array([Y_N, Y_O, Y_NO, Y_N2, Y_O2])
    # return rho_inf,

if __name__=="__main__":

    print('hi')

    mix = init_mixture()

    print("=============")
    print(" Inf state ")
    print("=============")
    mix.setState(rho_inf * Y_vec, T_inf, 1)
    a_inf     = mix.frozenSoundSpeed()
    u_scale = a_inf
    p_inf     = mix.P()
    Ma_inf    = u_inf / a_inf
    mu_inf    = mix.viscosity()
    kappa_inf = mix.frozenThermalConductivity()
    rhoE_inf = rho_inf * u_scale**2
    rho_vec_string = ' '.join(map(str, rho_inf * Y_vec))
    print("U0dim = [rho_i_inf, (rhou)_inf, (rhov)_inf, rhoE] = [" + rho_vec_string + " " + str(u_inf * rho_inf) + " " + str(0.0) + " " + str(rhoe(mix, u_inf)) + "]")
    print("U0nondim = [rho_i_inf, (rhou)_inf, (rhov)_inf, rhoE] = [" + " ".join(map(str, Y_vec)) + " " + str(u_inf * rho_inf / (rho_inf * a_inf)) + " " + str(0.0) + " " + str(rhoe(mix, u_inf) / rhoE_inf) + "]")

    print("PRINT SCALING TERMS")
    print("Just to make sure once again: (T, P) = (" + str(mix.T()) + ", " + str(mix.P()) + ")")

    print("rho_ref = " + str(rho_inf))
    print("u_ref = " + str(u_scale))
    print("rhoE_ref = " + str(rhoE_inf))
    print("T_ref = " + str(T_inf))
    print("mu_ref = " + str(mu_inf))
    print("kappa_ref = " + str(kappa_inf))

    print("MISC Terms")
    print("T_ref from cp: " + str(u_scale**2 / mix.mixtureFrozenCpMass()))
    print("Re =  " + str(rho_inf * u_scale / mu_inf))
    print("Pr = " + str(mu_inf * mix.mixtureFrozenCpMass() / kappa_inf))




    # print("PRINT SCALING TERMS using inf domain")
    # print("Just to make sure once again: (T, P) = (" + str(mix.T()) + ", " + str(mix.P()) + ")")
    # rho_inf = mix.density()
    # u_inf = mix.frozenSoundSpeed()
    # print("rho_inf = " + str(rho_inf))
    # print("u_inf = " + str(u_inf))
    # print("rhoE_inf = " + str(rho_inf * u_inf * u_inf))

    # cp = mix.mixtureFrozenCpMass()
    # gam = mix.mixtureFrozenGamma()
    # print("cp = " + str(cp))
    # print("gam = " + str(gam))

    # opt2 = mpp.MixtureOptions("air_5")
    # opt2.setStateModel("ChemNonEq1T")
    # opt2.setThermodynamicDatabase("NASA-9")
    # opt2.setViscosityAlgorithm("Gupta-Yos")
    # mix2 = mpp.Mixture(opt2)

    # print("=============")
    # print(" LEFT STATE ")
    # print("=============")
    # print_equil_vals(mix2, TL, PL)

    # print("=============")
    # print(" RIGHT STATE ")
    # print("=============")
    # print_equil_vals(mix2, TR, PR)

    # print()
    # print("PRINT SCALING TERMS using right side of domain")
    # print("Just to make sure once again: (T, P) = (" + str(mix2.T()) + ", " + str(mix2.P()) + ")")
    # rho_inf = mix2.density()
    # u_inf = mix2.frozenSoundSpeed()
    # print("rho_inf = " + str(rho_inf))
    # print("u_inf = " + str(u_inf))
    # print("rhoE_inf = " + str(rho_inf * u_inf * u_inf))

    # cp = mix2.mixtureFrozenCpMass()
    # gam = mix2.mixtureFrozenGamma()
    # print("cp = " + str(cp))
    # print("gam = " + str(gam))


