function out = mixtureFrozenCvMole(T, Mw, Y, X, thermo_structs)
    out =  mixtureFrozenCvMass(T, Mw, Y, thermo_structs) * mixtureMw(Mw, X);
end
