
T_dim = 1000*rand;
rho_i_dim = 1e-2*rand(5,1);
rhoe_dim = 1e-1*rand;

[dT_drhos1, dT_drhoe1, Dvec1, hvec1, mud1, kappa1, lambda1] = transportcoefficients(T_dim, rho_i_dim);
[dT_drhos2, dT_drhoe2, Dvec2, hvec2, mud2, kappa2, lambda2] = transportcoefficients2(T_dim, rho_i_dim);

max(abs(dT_drhos1-dT_drhos2))
max(abs(dT_drhoe1-dT_drhoe2))
max(abs(Dvec1-Dvec2))
max(abs(hvec1-hvec2))
max(abs(mud1-mud2))
max(abs(kappa1-kappa2))
max(abs(lambda1-lambda2))

kinetics_params = kinetics;
[species_thermo_structs, Mw, RU] = thermodynamicsModels;

omega1 = kineticsource(T_dim, rho_i_dim);
omega2 = netProductionRatesTotal(rho_i_dim, T_dim, Mw, kinetics_params, species_thermo_structs);    
max(abs(omega1-omega2))

f1 = equationofstate(T_dim, rho_i_dim, rhoe_dim);
rho_tilde = rho_i_dim ./ Mw;
alpha = -sum(rho_tilde);
f2 = f_T(T_dim, rho_tilde, rhoe_dim, alpha, species_thermo_structs);
max(abs(f1-f2))

