function out = viscosity(T, X_i, Mw, species_structs)
    mu_i = speciesViscosities(T, species_structs);
    phi_i = euckenPhi(mu_i, Mw, X_i);
    out = viscosityWilke(mu_i, X_i, phi_i);
end
