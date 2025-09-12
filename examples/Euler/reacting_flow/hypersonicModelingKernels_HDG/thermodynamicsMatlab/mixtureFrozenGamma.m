function out = mixtureFrozenGamma(T, Mw, Y, X, thermo_structs)
    cp = mixtureFrozenCpMole(T, Mw, Y, X, thermo_structs);
    cv = mixtureFrozenCvMole(T, Mw, Y, X, thermo_structs);
    out = cp / cv;
end
