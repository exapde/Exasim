function h_m = getEnthalpiesMass(T, Mw, species_thermo_structs)
    RU = 8.314471468617452;
    h = h_n(T, species_thermo_structs);
    h_m = reshape(h, size(Mw)) * T * RU ./ Mw;
end
