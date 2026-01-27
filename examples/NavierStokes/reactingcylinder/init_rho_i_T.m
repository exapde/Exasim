function [constants] = init_rho_i_T(rho_i, T, v_inf, L_ref)
    constants = struct();
    [species_thermo_structs, Mw, ~] = thermodynamicsModels();
    [blottner_structs, gupta_structs, gupta_mu_structs, gupta_kappa_structs] = transport();


    rho = sum(rho_i);
    Y = Y_i(rho_i);
    X = X_i(rho_i, Mw);
    P = pressure(T, rho_i, Mw);

    a = frozenSoundSpeed(T, rho_i, Mw, Y, X, species_thermo_structs);
    e = mixtureEnergyMass(T, P, X, Mw, species_thermo_structs);
    re = rho*e;
    rhoE = re + 0.5 * rho * v_inf^2;

    gam = mixtureFrozenGamma(T, Mw, Y, X, species_thermo_structs);
    cp = mixtureFrozenCpMass(T, Mw, Y, species_thermo_structs);

    mu_i = speciesViscosities(T, gupta_mu_structs);
    phi_i = euckenPhi(mu_i, Mw, X);
    mu_d = wilkeMixture(mu_i, X, phi_i);
    lambda_i = mu_d * 3/2 .* getCpsMass(T, Mw, species_thermo_structs);
    kappa = wilkeMixture(lambda_i, X, phi_i);

    
    constants.rho_ref = rho;
    constants.u_ref = a;    
    constants.rhoE_ref = rho * a^2; 
    constants.T_ref = T * (gam-1);
    constants.mu_ref = mu_d;
    constants.kappa_ref = kappa;
    constants.cp = cp; 
    constants.L_ref = L_ref;

    constants.rho_inf = rho_i;
    constants.rhou_inf = rho * v_inf;
    constants.rhov_inf = 0.0;
    constants.rhoE_inf = rhoE;

end