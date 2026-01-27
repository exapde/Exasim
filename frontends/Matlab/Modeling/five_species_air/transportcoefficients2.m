function [dT_drho_i_dim, dT_drhoe_dim, D_vec, h_vec, mu_d_dim, kappa_dim, lambda_i] = transportcoefficients2(T_dim, rho_i_dim)

    [species_thermo_structs, Mw, RU] = thermodynamicsModels();
    [~, gupta_structs, gupta_mu_structs, ~] = transport();
    
    Y = Y_i(rho_i_dim);
    denom = sum(rho_i_dim) * mixtureFrozenCvMass(T_dim, Mw, Y, species_thermo_structs);
    e_i = getEnergiesMass(T_dim, Mw, species_thermo_structs);
    dT_drho_i_dim = -e_i(:) ./ denom;
    dT_drhoe_dim = 1.0 / denom;

    p_dim = pressure(T_dim, rho_i_dim, Mw);
    X = X_i(rho_i_dim,Mw);
    D_vec = averageDiffusionCoeffs(T_dim, X, Y, Mw, p_dim, gupta_structs);
    h_vec = getEnthalpiesMass(T_dim, Mw, species_thermo_structs);
      
    mu_i = speciesViscosities(T_dim, gupta_mu_structs);
    phi_i = euckenPhi(mu_i, Mw, X);
    mu_d_dim = wilkeMixture(mu_i, X, phi_i);

    lambda_i = mu_d_dim * (getCpsMass(T_dim, Mw, species_thermo_structs) + 5/4 * RU./Mw);
    kappa_dim = wilkeMixture(lambda_i, X, phi_i);
end
