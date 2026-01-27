 function [rho, rhoe_dim] = density_energy_from_YiTP(Y_i, T, p)

    [species_thermo_structs, Mw, RU] = thermodynamicsModels();
    X = (Y_i .* Mw) / sum(Y_i .* Mw);
    rho = density(T, p, X, Mw);
    rho_i = rho * Y_i;
    e_dim = mixtureEnergyMass(T, p, X, Mw, species_thermo_structs);
    rhoe_dim = (sum(rho) * e_dim);
end