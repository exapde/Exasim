function out = mixtureEnergyMass(T, P, X, Mw, species_thermo_structs)
    RU = 8.314471468617452;
    rho = density(T, P, X, Mw);
    hmass = mixtureHMass(T, X, Mw, species_thermo_structs);
    out =  hmass - P ./ rho;
end