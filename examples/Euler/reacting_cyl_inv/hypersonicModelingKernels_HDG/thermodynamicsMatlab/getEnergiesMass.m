function e_i = getEnergiesMass(T,Mw, species_thermo_structs)
    RU = 8.314471468617452;
    h = h_n(T, species_thermo_structs);
    if size(h,2) ~= size(Mw,2)
        h = h';
    end
    e_i = (h - 1.0)*T*RU ./ Mw;
end
