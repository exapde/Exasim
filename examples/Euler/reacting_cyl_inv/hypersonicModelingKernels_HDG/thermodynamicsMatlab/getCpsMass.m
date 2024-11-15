function cp_m = getCpsMass(T, Mw, species_thermo_structs)
    RU = 8.314471468617452;
    cp = cp_n(T, species_thermo_structs);
    cp_m = reshape(cp, size(Mw)) * RU ./ Mw;
end
