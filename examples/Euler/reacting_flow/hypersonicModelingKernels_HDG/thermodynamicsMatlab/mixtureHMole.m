function out = mixtureHMole(T, X, species_thermo_structs)
    RU = 8.314471468617452;
    h = h_n(T, species_thermo_structs);
    hmole = sum( h .* X);
    hmole = hmole * RU * T;
    out =  hmole;
end
