function out = mixtureFrozenCpMass(T, Mw, Y, species_thermo_structs)
    cp = getCpsMass(T, Mw, species_thermo_structs);
    out = sum(cp .* Y);
end
