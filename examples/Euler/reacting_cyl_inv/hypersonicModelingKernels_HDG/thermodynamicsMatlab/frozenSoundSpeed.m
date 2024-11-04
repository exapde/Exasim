function out = frozenSoundSpeed(T, r_i, Mw, Y, X, thermo_structs)
    RU = 8.314471468617452;
    gam = mixtureFrozenGamma(T, Mw, Y, X, thermo_structs);
    P = pressure(T, r_i, Mw);
    density = sum(r_i);
    out = sqrt(gam * P / density);
end