function cv_m = getCvsMass(T, Mw, species_thermo_structs)
    RU = 8.314471468617452;
    cp = cp_n(T, species_thermo_structs);
    if size(cp,2) ~= size(Mw,2)
        cp = cp';
    end
    cv_m = (cp - 1.0) * RU ./ Mw;
end
