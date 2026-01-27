function out = mixtureHMass(T, X, Mw, species_thermo_structs)
    hmole = mixtureHMole(T, X, species_thermo_structs);
    out =  hmole ./ mixtureMw(Mw, X);
end
