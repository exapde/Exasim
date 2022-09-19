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

    print("To check")
    print("(T, P) = (" + str(mix.T()) + ", " + str(mix.P()) + ")")
    print("Initial conditions (dimensional): ")
    print("(rho_i, rhou, rhoE) = (" + str(rhovec) + " 0.0 " + str(rhoE) + ")")

def init_mixture():
    opt = mpp.MixtureOptions("air_5")
    opt.setStateModel("ChemNonEq1T")
    opt.setViscosityAlgorithm("Gupta-Yos")

    mix = mpp.Mixture(opt)
    return mix

def test_range():
    TL = 300.0
    PL = 1.0e4

    opt = mpp.MixtureOptions("air_5")
    opt.setStateModel("ChemNonEq1T")
    opt.setViscosityAlgorithm("Gupta-Yos")
    opt.setThermodynamicDatabase("RRHO")

    mix = mpp.Mixture(opt)

    for i in range(0,1000):
        TL += i
        mix.equilibrate(TL, PL)
        print(mix.mixtureEnergyMass())




if __name__=="__main__":
    TL = 9000.0 + 298.17
    PL = 1.95256e5

    TR = 300.0 + 298.17
    PR = 1.0e4

    opt = mpp.MixtureOptions("air_5")
    opt.setStateModel("ChemNonEq1T")
    opt.setViscosityAlgorithm("Gupta-Yos")
    opt.setTeh

    mix = mpp.Mixture(opt)

    print("=============")
    print(" LEFT STATE ")
    print("=============")
    print_equil_vals(mix, TL, PL)

    print("=============")
    print(" RIGHT STATE ")
    print("=============")
    print_equil_vals(mix, TR, PR)

    print()
    print("PRINT SCALING TERMS using right side of domain")
    print("Just to make sure once again: (T, P) = (" + str(mix.T()) + ", " + str(mix.P()) + ")")
    rho_inf = mix.density()
    u_inf = mix.frozenSoundSpeed()
    print("rho_inf = " + str(rho_inf))
    print("u_inf = " + str(u_inf))
    print("rhoE_inf = " + str(rho_inf * u_inf * u_inf))

    cp = mix.mixtureFrozenCpMass()
    gam = mix.mixtureFrozenGamma()
    print("cp = " + str(cp))
    print("gam = " + str(gam))


