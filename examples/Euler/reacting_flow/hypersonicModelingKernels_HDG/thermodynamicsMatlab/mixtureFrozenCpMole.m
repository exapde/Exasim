function out = mixtureFrozenCpMole(T, Mw, Y, X, thermo_structs)
    out =  mixtureFrozenCpMass(T, Mw, Y, thermo_structs) * mixtureMw(Mw, X);
end
